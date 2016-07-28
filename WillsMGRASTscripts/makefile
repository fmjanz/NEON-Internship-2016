
TARGET=$(shell cat targets2)
RSPHYLUM=$(addsuffix .RS.phylum.csv,$(TARGET) ) 
RDPPHYLUM=$(addsuffix .RDP.phylum.csv,$(TARGET)) 
RDPSTRAIN= $(addsuffix .RDP.strain.csv,$(TARGET))
RSSTRAIN=$(addsuffix .RS.strain.csv,$(TARGET)) 
SSFUNCTION=$(addsuffix .SS.function.csv,$(TARGET))  
JSON=$(addsuffix .json,$(TARGET))  
ARCHAIC= $(addsuffix .refseq.level6.csv,$(TARGET)) $(addsuffix .ss.level3.csv,$(TARGET))
ARCH2= $(addsuffix .RS.species.csv,$(TARGET)) $(addsuffix .SS.function.csv,$(TARGET))

rdpphylum: $(RDPPHYLUM)
rsphylum: $(RSPHYLUM)
rdpstrain: $(RDPSTRAIN)
rsstrain: $(RSSTRAIN)
ssfunction: $(SSFUNCTION)
overview: $(JSON)
archaic: $(ARCHAIC)
archaic2: $(ARCH2)

all: $(JSON) $(RDPSTRAIN) $(RSSTRAIN) $(SSFUNCTION)

%.RS.phylum.csv: 
	-{ time ./abundance-matrix.py -m $* -s RefSeq -g phylum   > $@ 2> $@.err; }  2> $@.time

%.RS.strain.csv: 
	-{ time ./abundance-matrix.py -m $* -s RefSeq -g strain   > $@ 2> $@.err; }  2> $@.time

%.RDP.strain.csv: 
	-{ time ./abundance-matrix.py -m $* -s RDP -g strain   > $@ 2> $@.err; }  2> $@.time

%.RDP.phylum.csv: 
	-{ time ./abundance-matrix.py -m $* -s RDP -g phylum   > $@ 2> $@.err; }  2> $@.time

%.json: 
	-{ curl 'http://api.metagenomics.anl.gov/1/metagenome/$*?verbosity=full' > $@; }

%.json.profile.csv: %.json
	-./parse-profiles.py $^ > $@

%.json.subsystem.csv: %.json
	-./parse-annotations.py $^ > $@

%.json.all.csv: %.json
	-./parse-all.py  $^ > $@

%.SS.function.csv: 
	-{ time ./abundance-matrix.py -m $* -c function -s Subsystems -g function > $@ 2> $@.err; }  2> $@.time

%.RS.species.csv: 
	-{ time ./abundance-matrix.py -m $* -c organism -s RefSeq -g species > $@ 2> $@.err; }  2> $@.time

%.refseq.level6.csv:
	-{ time ./abundance-matrix.py RefSeq  species $* > $@ 2> $@.err ; }

%.ss.level3.csv:
	-{ time ./abundance-matrix.py Subsystems level3 $* > $@ 2> $@.err ; }

clean:
	ls | grep err  | xargs grep ERROR | grep -v "no data"  > ERRORS
	-mkdir broken
	mv $(shell cat ERRORS | sed 's/.err/* /' | cut -f 1 -d ' ' ) broken
#	mv $(shell cat ERRORS | sed 's/.err/.time /' | cut -f 1 -d ' ' ) broken
#	mv $(shell cat ERRORS | sed 's/.err/.err /' | cut -f 1 -d ' ' ) broken

