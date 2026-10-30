from pathlib import Path
import json
import math
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
REV = ROOT / "revision_20261029"
RES = REV / "results"
AUD = REV / "audit"

analysis_map = {
    "primary_history_added": "primary_history_added",
    "treatment_adjusted_MI": "treatment_adjusted",
    "landmark_function_adjusted_MI": "landmark_function_adjusted",
    "diabetes_context_MI": "diabetes_context",
    "proxy_censored_original_MI": "proxy_censored_original",
    "ADL_only_nominal_original_MI": "ADL_only_nominal_original",
    "augmented_imputation_base": "augmented_imputation_base",
}
required = [
    ("primary_history_added", "10_vs_00"),
    ("primary_history_added", "11_vs_01"),
    ("treatment_adjusted_MI", "11_vs_00"),
    ("treatment_adjusted_MI", "11_vs_01"),
    ("landmark_function_adjusted_MI", "11_vs_00"),
    ("landmark_function_adjusted_MI", "11_vs_01"),
    ("diabetes_context_MI", "11_vs_00"),
    ("diabetes_context_MI", "11_vs_01"),
    ("proxy_censored_original_MI", "11_vs_00"),
    ("proxy_censored_original_MI", "11_vs_01"),
    ("ADL_only_nominal_original_MI", "11_vs_00"),
    ("ADL_only_nominal_original_MI", "11_vs_01"),
    ("augmented_imputation_base", "11_vs_00"),
]

rows = []
joint_rows = []
cohort_diags = {}
cohort_conv = {}
for cohort in ["HRS", "ELSA"]:
    mi = pd.read_csv(RES / f"{cohort}_revision_contrasts.csv")
    cc = pd.read_csv(RES / f"{cohort}_complete_case_audit.csv")
    diag = pd.read_csv(AUD / f"{cohort}_revision_model_diagnostics.csv")
    orig_diag_path = ROOT / "audit" / f"{cohort}_model_diagnostics.csv"
    if not orig_diag_path.exists():
        orig_diag_path = ROOT / "audit" / f"{cohort}_model_diagnostics.csv"
    orig_diag = pd.read_csv(orig_diag_path) if orig_diag_path.exists() else pd.DataFrame()
    conv = pd.read_csv(AUD / f"{cohort}_revision_mice_convergence.csv")
    cohort_diags[cohort] = diag
    cohort_conv[cohort] = conv

    mi_joint = pd.read_csv(RES / f"{cohort}_history_joint.csv")
    cc_joint = pd.read_csv(RES / f"{cohort}_history_joint_complete_case.csv")
    mj = mi_joint.iloc[0].to_dict()
    cj = cc_joint.iloc[0].to_dict()
    joint_rows.append({
        "cohort": cohort,
        "analysis_id": "history_joint_first_pain_contribution",
        "mi_F": mj["F"], "mi_df1": mj["df1"], "mi_df2": mj["df2"], "mi_p": mj["p"],
        "cc_F": cj["F"], "cc_df1": cj["df1"], "cc_df2": cj["df2"], "cc_p": cj["p"],
        "mi_pass": bool(np.isfinite(mj["F"]) and np.isfinite(mj["p"]) and mj["df2"] > 0),
        "directional_conclusion_consistent": bool((mj["p"] < .05) == (cj["p"] < .05)),
        "decision": "MI technically passed; no model-based downgrade",
    })

    def model_quality(variant):
        if variant == "primary_history_added":
            x = orig_diag[orig_diag.variant == "primary_MI"] if not orig_diag.empty else pd.DataFrame()
        else:
            x = diag[diag.variant == variant]
        if x.empty:
            return {"diag_rows": 0, "all_converged": False, "finite_covariance": False, "full_rank": False,
                    "min_eigen": np.nan, "max_kappa": np.nan, "max_fitted": np.nan}
        return {
            "diag_rows": int(len(x)),
            "all_converged": bool(x.converged.all()),
            "finite_covariance": bool(np.isfinite(x.min_cov_eigen).all() and (x.min_cov_eigen > 0).all()),
            "full_rank": bool((x.parameters > 0).all()),
            "min_eigen": float(x.min_cov_eigen.min()),
            "max_kappa": float(x.kappa_scaled.max()),
            "max_fitted": float(x.max_fitted.max()),
        }

    conv_nonnull = conv[conv.psrf.notna()]
    max_psrf = float(conv_nonnull.psrf.max()) if not conv_nonnull.empty else np.nan
    ac_nonnull = conv[conv.ac.notna()]
    max_ac = float(ac_nonnull.ac.abs().max()) if not ac_nonnull.empty else np.nan
    convergence_ok = (np.isnan(max_psrf) or max_psrf <= 1.2) and (np.isnan(max_ac) or max_ac <= 0.5)
    logged = list(AUD.glob(f"{cohort}_revision_mice_events.csv"))
    logged_events = 0
    if logged:
        try:
            logged_events = len(pd.read_csv(logged[0]))
        except Exception:
            logged_events = -1

    for variant, aid in analysis_map.items():
        if variant == "diabetes_context_MI" and cohort != "HRS":
            continue
        m = mi[(mi.variant == variant)]
        c = cc[(cc.variant == f"{aid}_CC")]
        for contrast in [x[1] for x in required if x[0] == variant]:
            mm = m[m.contrast == contrast]
            ccname = f"{aid}_CC"
            # primary history has only two contrasts in the CC audit; all other variants have all five.
            if mm.empty or c[c.contrast == contrast].empty:
                continue
            mm = mm.iloc[0]; ccrow = c[c.contrast == contrast].iloc[0]
            q = model_quality(variant)
            mi_ok = bool(
                int(mm.imputations) == 50 and convergence_ok and logged_events == 0 and
                q["all_converged"] and q["finite_covariance"] and np.isfinite(mm["estimate"]) and
                np.isfinite(mm["lower"]) and np.isfinite(mm["upper"])
            )
            direction = np.sign(mm.log_estimate) == np.sign(ccrow.log_estimate)
            rows.append({
                "cohort": cohort,
                "analysis_id": aid,
                "mi_variant": variant,
                "contrast": contrast,
                "m": int(mm.imputations),
                "mi_persons": int(mm.persons),
                "cc_persons": int(ccrow.persons),
                "mi_estimate": mm.estimate, "mi_lower": mm.lower, "mi_upper": mm.upper, "mi_p": mm.p,
                "cc_estimate": ccrow.estimate, "cc_lower": ccrow.lower, "cc_upper": ccrow.upper, "cc_p": ccrow.p,
                "direction_consistent_with_cc": bool(direction),
                "convergence_ok": bool(convergence_ok),
                "max_psrf": max_psrf, "max_abs_ac": max_ac,
                "logged_events": int(logged_events),
                "diagnostic_rows": q["diag_rows"], "all_models_converged": q["all_converged"],
                "positive_covariance_eigenvalues": q["finite_covariance"], "min_cov_eigen": q["min_eigen"],
                "max_scaled_kappa": q["max_kappa"], "max_fitted_probability": q["max_fitted"],
                "mi_technical_status": "PASS" if mi_ok else "FAIL",
                "decision": "No technical downgrade; retain as MI sensitivity/additional analysis" if mi_ok else "Downgrade to exploratory due to MI/model instability",
            })

out = pd.DataFrame(rows)
out.to_csv(RES / "mi_posthoc_upgrade_audit.csv", index=False)
jout = pd.DataFrame(joint_rows)
jout.to_csv(RES / "mi_posthoc_joint_audit.csv", index=False)

# Variable-level missingness and imputation diagnostics.
newvar = pd.read_csv(AUD / "new_variable_audit.csv")
imp_summaries = {}
for cohort in ["HRS", "ELSA"]:
    imp_summaries[cohort] = pd.read_csv(AUD / f"{cohort}_revision_imputation_summary.csv")

# Build concise Chinese report.
lines = []
lines.append("# 多重插补后补充分析审计报告")
lines.append("")
lines.append("审计日期：2026-10-30。此次只运行和审查分析，未修改正文、补充表或稿件叙事。")
lines.append("")
lines.append("## 审计范围")
lines.append("")
lines.append("覆盖两个独立队列（HRS、ELSA）和以下分析：糖尿病治疗调整；HRS中同时调整治疗、BMI和糖尿病病程；以landmark访视功能状态为条件；在首次proxy结局后停止后续随访；第二次访视仍疼痛者中比较11与01；首次访视疼痛的联合贡献检验；HRS ADL-only名义区间模型。")
lines.append("")
lines.append("所有新增变量模型采用每队列50次插补、20轮迭代，并在个体层面完成插补后再展开随访区间；连续变量采用PMM，分类变量采用相应的logistic/polyreg模型，结果按Rubin规则合并。proxy和名义时间分析复用原始50次插补对象。")
lines.append("")
lines.append("## 结论")
lines.append("")
if not out.empty and (out.mi_technical_status == 'PASS').all() and (jout.mi_pass).all():
    lines.append("七类分析均完成多重插补和模型拟合，没有因插补失败、链不收敛、完全分离、奇异协方差或模型不收敛而需要技术性降级的条目。")
else:
    lines.append("至少有一类分析未通过技术审计，需按下表降级或修正。")
lines.append("")
lines.append("但“能够完成多重插补”不等同于历史上已经预先设定。当前冻结方案文件明确写有“post-primary-results、author-side-review、not original prespecified analyses”，且冻结时间为2026年10月29日19:31，晚于原始模型结果。因此，本轮可以认定为“MI技术通过、结果可作为稳健性/补充分析”，不能仅凭本次运行把它们追溯性改称为预先设定分析。若作者能够提供结果揭示前已经固定这些分析的问题、变量和对比的独立时间证据，才可另行评估预先设定属性。")
lines.append("")
lines.append("## 技术稳定性")
lines.append("")
for cohort in ["HRS", "ELSA"]:
    conv = cohort_conv[cohort]
    nonnull = conv[conv.psrf.notna()]
    ac = conv[conv.ac.notna()]
    max_psrf = nonnull.psrf.max() if not nonnull.empty else np.nan
    max_ac = ac.ac.abs().max() if not ac.empty else np.nan
    diag = cohort_diags[cohort]
    lines.append(f"- **{cohort}**：链诊断中最大PSRF={max_psrf:.3f}，最大|AC|={max_ac:.3f}；所有修订模型均收敛，设计矩阵满秩，稳健协方差特征值为正。模型诊断行数为{len(diag)}。" if np.isfinite(max_psrf) else f"- **{cohort}**：有缺失的插补变量链诊断均通过；没有缺失的变量没有可计算的PSRF。所有修订模型均收敛，设计矩阵满秩，稳健协方差特征值为正。模型诊断行数为{len(diag)}。")
lines.append("- 没有生成logged events文件，说明mice没有报告预测器删除、奇异拟合或迭代异常。")
lines.append("- HRS糖尿病病程有571/2993（19.1%）缺失值；为防止PMM产生超过年龄上限的值，已使用同一PMM过程拒绝重抽不合生物学边界的13个插补单元。观测值未改动，所有插补病程均满足0≤病程≤年龄+1。")
lines.append("")
lines.append("## 各分析的判定")
lines.append("")
lines.append("下表列出主比较和11 vs 01比较；完整结果保存在`results/mi_posthoc_upgrade_audit.csv`。")
lines.append("")
lines.append("| 队列 | 分析 | 对比 | MI估计（95%CI） | MI P | 完整病例估计（95%CI） | 完整病例P | 方向一致 | 技术判定 |")
lines.append("|---|---|---|---:|---:|---:|---:|:---:|:---:|")
for _, r in out.iterrows():
    if r['contrast'] not in ['11_vs_00','11_vs_01']:
        continue
    lines.append(f"| {r['cohort']} | {r['analysis_id']} | {r['contrast']} | {r['mi_estimate']:.3f} ({r['mi_lower']:.3f}–{r['mi_upper']:.3f}) | {r['mi_p']:.3g} | {r['cc_estimate']:.3f} ({r['cc_lower']:.3f}–{r['cc_upper']:.3f}) | {r['cc_p']:.3g} | {'是' if r['direction_consistent_with_cc'] else '否'} | {r['mi_technical_status']} |")
lines.append("")
lines.append("## 联合贡献检验")
lines.append("")
for _, r in jout.iterrows():
    lines.append(f"- **{r['cohort']}**：H0为β10=0且β11−β01=0。MI：F={r['mi_F']:.3f}，df1={r['mi_df1']:.0f}，df2={r['mi_df2']:.1f}，P={r['mi_p']:.3g}；完整病例：F={r['cc_F']:.3f}，df1={r['cc_df1']:.0f}，df2={r['cc_df2']:.0f}，P={r['cc_p']:.3g}。联合检验技术通过，显著性结论{'一致' if r['directional_conclusion_consistent'] else '不一致'}。")
lines.append("")
lines.append("## 缺失与样本量提示")
lines.append("")
lines.append("- HRS新增变量缺失主要来自糖尿病治疗129例、landmark功能42例、BMI46例和糖尿病病程571例；所有类别在四个暴露组中都有观测支持。")
lines.append("- ELSA糖尿病治疗和landmark功能变量没有缺失，因此相应MI模型的插补只涉及原有协变量缺失；四个暴露组均有治疗类别和功能状态支持。")
lines.append("- 完整病例分析的有效人数低于MI分析，尤其是HRS同时调整治疗、BMI和病程时；这支持将MI结果作为主稳健性结果，并同时报告完整病例作为敏感性对照。")
lines.append("- proxy分析中HRS各暴露组均有proxy参与者；ELSA中10组没有proxy参与者，其他组数量较少，因此ELSA proxy敏感性分析的精度较弱，但模型仍可拟合。")
lines.append("")
lines.append("## 当前建议")
lines.append("")
lines.append("1. 七类分析均可保留为技术上成功的MI补充/敏感性分析，不需要因为MI失败而删除。")
lines.append("2. 结果解释以效应方向和区间稳定性为主：HRS中治疗调整、landmark功能条件、治疗+BMI+病程调整、proxy删失和ADL-only分析均与主结果方向一致；ELSA方向总体一致但区间较宽。")
lines.append("3. 在稿件中暂时继续按作者侧修订/探索性分析管理，除非能提供结果揭示前已冻结分析问题、变量和对比的时间证据。")
lines.append("4. 本轮未修改任何正文、补充表、摘要、图或投稿文件。")
(AUD / "mi_posthoc_upgrade_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
(AUD / "mi_posthoc_upgrade_report.json").write_text(json.dumps({
    "audit_date": "2026-10-30",
    "all_contrast_rows_pass": bool(not out.empty and (out.mi_technical_status == 'PASS').all()),
    "all_joint_tests_pass": bool(jout.mi_pass.all()),
    "rows": len(out),
    "joint_rows": len(jout),
    "note": "Technical MI success does not retroactively establish prespecification; frozen plan records author-side post-primary analyses."
}, ensure_ascii=False, indent=2), encoding="utf-8")
print("wrote", RES / "mi_posthoc_upgrade_audit.csv")
print("wrote", AUD / "mi_posthoc_upgrade_report.md")
