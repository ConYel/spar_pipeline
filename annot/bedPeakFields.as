table hg19sncRNA
"DASHR sncRNA annotation"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Peak"
uint    score;		"Score"
char[1] strand;		"Strand"
string  rnaID;	"DASHR annotation"
float  readCount;	"Expression [raw read count]"
uint  maxPos5p; 	"Position with most 5p read ends"
uint  maxPos3p; 	"Position with most 3p read ends"
float  maxProb5p; 	"Fraction of reads with same 5p end"
float maxProb3p; 	"Fraction of reads with same 3p end"
)
