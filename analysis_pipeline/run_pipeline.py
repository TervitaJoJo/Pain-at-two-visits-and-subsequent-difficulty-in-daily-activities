"""Portable analysis runner. All generated participant data stay in --work-dir."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys,time,hashlib
HERE=Path(__file__).resolve().parent
CODE=HERE.parent
OLD=Path('analysis/pain_transition_diabetes_adl_20261029')
EXT=Path('analysis/pain_external_replication_20261030')
ADD=EXT/'additional_20261030';REV=OLD/'revision_20261029';META=EXT/'meta_analysis'
def write(p,text):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text,encoding='utf-8',newline='\n')
def stage_sources(work):
 for src in (HERE/'source').rglob('*'):
  if src.is_file():
   dst=work/src.relative_to(HERE/'source');dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
 for base in (OLD,REV,EXT,ADD,META):
  for folder in ('data_private','results','audit','reports','figures'):(work/base/folder).mkdir(parents=True,exist_ok=True)
 for folder in ('qa','main_figures','supplementary_figures'):(work/'presentation'/folder).mkdir(parents=True,exist_ok=True)
def metadata(data,work):
 import pyreadstat
 paths={'HRS':data/'HRS/RAND HRS Data/Longitudinal and Cross-Wave Data Products/randhrs1992_2022v1.dta',
        'ELSA':data/'ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/h_elsa_g3.dta'}
 records={}
 for cohort,p in paths.items():
  if not p.is_file():raise FileNotFoundError(f'Required authorized input missing: {p}')
  _,m=pyreadstat.read_dta(p,metadataonly=True)
  records[cohort]={'path':str(p),'labels':m.column_names_to_labels,'columns':[{'name':n,'label':m.column_names_to_labels.get(n,'')} for n in m.column_names]}
 write(work/'analysis/second_line_rescreen_20261028/local/community_metadata.json',json.dumps({'ELSA_harmonized':records['ELSA']},indent=2))
 write(work/'analysis/cross_system_multimorbidity_mortality_20261028/results/metadata_HRS_rand.json',json.dumps(records['HRS'],indent=2))
def run(cmd,work,env):
 stamp=time.time();name=f'{len(list((work/"execution_logs").glob("*.json")))+1:03d}_{Path(cmd[1]).stem}'
 print('RUN', ' '.join(map(str,cmd)),flush=True)
 result=subprocess.run(list(map(str,cmd)),cwd=work,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
 write(work/'execution_logs'/f'{name}.log',result.stdout+'\n'+result.stderr)
 write(work/'execution_logs'/f'{name}.json',json.dumps({'command':list(map(str,cmd)),'returncode':result.returncode,'seconds':time.time()-stamp},indent=2))
 if result.returncode:raise RuntimeError(f'Exit {result.returncode}; see execution_logs/{name}.log')
 if result.stdout.strip():print(result.stdout[-1500:],flush=True)
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--work-dir',type=Path,required=True,help='External writable run directory; must not be inside this package.')
 p.add_argument('--data-root',type=Path,help='Authorized raw cohort data root, required for prepare/all.')
 p.add_argument('--rscript',default='Rscript')
 p.add_argument('--stage',choices=['all','prepare','primary','revision','external','additional','synthesis','verify','figures','summary-check'],default='all')
 a=p.parse_args();work=a.work_dir.resolve()
 if work==CODE or CODE in work.parents or work in CODE.parents:raise ValueError('work-dir must be separate from the frozen code tree')
 marker=work/'.anonymous_analysis_workspace'
 if work.exists() and any(work.iterdir()) and not marker.exists():raise ValueError('Choose an empty work-dir or a workspace created by this runner')
 work.mkdir(parents=True,exist_ok=True);marker.write_text('Generated analysis workspace; may contain restricted participant data.\n',encoding='utf-8');(work/'execution_logs').mkdir(exist_ok=True)
 stage_sources(work)
 env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONIOENCODING='utf-8',MPLBACKEND='Agg',REUSE_MI_CACHE='0')
 env['PYTHONPATH']=str(work/'presentation')
 if a.data_root:env['COHORT_DATA_ROOT']=str(a.data_root.resolve())
 py=lambda script,*args:run([sys.executable,str(script),*args],work,env)
 rr=lambda script,*args:run([a.rscript,str(script),*args],work,env)
 if a.stage=='summary-check':
  # No raw data or fitted-model caches are required. References are copied only here.
  for rec in json.loads((HERE/'reference_index.json').read_text()):
   dst=work/rec['stage_to']
   if dst.exists():raise FileExistsError('summary-check requires a fresh workspace to avoid replacing computed results: '+str(dst))
   dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(CODE/rec['file'],dst)
  rr(META/'code/01_meta.R',META);py(META/'code/02_verify.py')
  rr(ADD/'code/05_meta_additional.R');py(ADD/'code/07_verify.py')
  py(ADD/'figures/plot_absolute_risks.py')
  py(META/'figures/plot_meta_forest.py')
  for f in sorted((work/'presentation/figures').glob('gen_*.py')):py(f.relative_to(work))
  print('Summary-only verification complete; private-data covariance checks were not run.');return
 stages=['prepare','primary','revision','external','additional','synthesis','verify','figures'] if a.stage=='all' else [a.stage]
 for stage in stages:
  if stage=='prepare':
   if a.data_root is None:p.error('--data-root is required for prepare/all')
   metadata(a.data_root.resolve(),work)
   py(OLD/'code/01_build_cohorts.py');py(REV/'01_extract_revision.py')
   py(REV/'audit_variables.py');py(REV/'04_attrition_by_pattern.py')
   py(EXT/'code/01_build.py');py(ADD/'code/01_extract_share.py')
  elif stage=='primary':
   for cohort in ('HRS','ELSA'):rr(OLD/'code/02_analyze.R',cohort)
   rr(OLD/'code/03_export_verification.R');py(OLD/'code/04_independent_verification.py');py(OLD/'code/05_compile_tables.py')
  elif stage=='revision':
   for cohort in ('HRS','ELSA'):rr(REV/'02_analyze_revision.R',cohort);rr(REV/'06_complete_case_audit.R',cohort)
   rr(REV/'03_verify_revision.R');py(REV/'07_mi_posthoc_upgrade_audit.py')
  elif stage=='external':
   for cohort in ('SHARE','CHARLS'):rr(EXT/'code/02_analyze.R',cohort)
   py(EXT/'code/03_verify_models.py');py(EXT/'code/04_compile.py')
  elif stage=='additional':
   rr(ADD/'code/02_share_additional.R');rr(ADD/'code/03_standardized_risks.R');rr(ADD/'code/06_export_verification.R')
  elif stage=='synthesis':rr(META/'code/01_meta.R',META);rr(ADD/'code/05_meta_additional.R')
  elif stage=='verify':py(META/'code/02_verify.py');py(ADD/'code/07_verify.py')
  elif stage=='figures':
   py(ADD/'figures/plot_absolute_risks.py');py(META/'figures/plot_meta_forest.py')
   for f in sorted((work/'presentation/figures').glob('gen_*.py')):py(f.relative_to(work))
 print('Completed:',','.join(stages))
if __name__=='__main__':main()
