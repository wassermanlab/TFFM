import sys
sys.path.append("/home/amathelier/PostDoc/HMM/GHMM/src/TFFM-framework_github")
import tffm_module
from constants import TFFM_KIND

tffm_first_order = tffm_module.tffm_from_meme("meme.txt", TFFM_KIND.FIRST_ORDER)
tffm_first_order.write("tffm_first_order_initial.xml")
tffm_first_order.train("train.fa")
tffm_first_order.write("tffm_first_order.xml")
out = open("tffm_first_order_summary_logo.svg", "w")
tffm_first_order.print_summary_logo(out)
out.close()
out = open("tffm_first_order_dense_logo.svg", "w")
tffm_first_order.print_dense_logo(out)
out.close()

tffm_detailed = tffm_module.tffm_from_meme("meme.txt", TFFM_KIND.DETAILED)
tffm_detailed.write("tffm_detailed_initial.xml")
tffm_detailed.train("train.fa")
tffm_detailed.write("tffm_detailed.xml")
out = open("tffm_detailed_summary_logo.svg", "w")
tffm_detailed.print_summary_logo(out)
out.close()
out = open("tffm_detailed_dense_logo.svg", "w")
tffm_detailed.print_dense_logo(out)
out.close()

tffm_first_order = tffm_module.tffm_from_xml("tffm_first_order.xml",
        TFFM_KIND.FIRST_ORDER)
print "1st-order all"
for hit in tffm_first_order.scan_sequences("test.fa"):
    if hit:
        print hit

print "1st-order best"
for hit in tffm_first_order.scan_sequences("test.fa", only_best=True):
    print hit

tffm_detailed = tffm_module.tffm_from_xml("tffm_detailed.xml",
        TFFM_KIND.DETAILED)
print "detailed all"
for hit in tffm_detailed.scan_sequences("test.fa"):
    if hit:
        print hit

print "detailed best"
for hit in tffm_detailed.scan_sequences("test.fa", only_best=True):
    print hit
