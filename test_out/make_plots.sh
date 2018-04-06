set -xv
/appl/R-3.2.3/bin/Rscript test_out/module3_annot.r >> test_out/logs/R.log
bash ./scripts/make_submodule_html.sh ./scripts/R/module3 test_out/figures > test_out/module3_annot.html
