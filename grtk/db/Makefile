DB_HOST=wrendb
DB_SCHEMA=gfa

psql=psql -h $(DB_HOST) -d $(DB_SCHEMA) -t --no-align -F '	'
NCBI_FTP=ftp://ftp.ncbi.nlm.nih.gov
DATA=data/gene_ontology_ext.obo \
	 data/taxdump.tar.gz \
	 data/gene2go.gz \
	 data/gene_info.gz

$(shell mkdir -p data/)

all: load-base load-go

load-base: $(DATA)
	$(psql) < sql/schema.sql

	tar Oxzf data/taxdump.tar.gz names.dmp | grep scientific \
		| cut -f1,3 \
		| $(psql) -c "COPY taxon FROM STDIN DELIMITER '	' NULL '-'"

	zcat data/gene_info.gz | grep -v "^#" | cut -f1-3,12 \
		| awk '$$1 != 366646' \
		| $(psql) -c "COPY gene FROM STDIN DELIMITER '	' NULL '-'" 


load-go: data/gene_ontology_ext.obo
	./script/import-obo $(DB_SCHEMA) $(DB_HOST) \
		GO "Gene Ontology" data/gene_ontology_ext.obo

	$(psql) -c "SELECT accession, term.id FROM term \
		INNER JOIN ontology \
		ON term.ontology_id=ontology.id \
		WHERE ontology.namespace='GO'" \
		| sort -k1,1 > data/goacc2termid

	zcat data/gene2go.gz | cut -f2-4 | sort -k2,2 \
		| join -1 2 -2 1 -o "1.1 2.2 1.3" - data/goacc2termid \
		| grep -v IEA | awk '{print $$1,$$2,"1"}' \
		| sort -u \
		| $(psql) -c "COPY gene_term FROM STDIN DELIMITER ' '"


load-brenda: data/manual_annotations_ursa.csv data/brenda.obo
	./script/import-obo $(DB_SCHEMA) $(DB_HOST) \
		BTO "BRENDA Tissue Ontology" data/brenda.obo
	./script/parse-ursa-annotations $< \
		| $(psql) -c \
		"CREATE TEMPORARY TABLE \
			ursa_annotations(gsm VARCHAR, tissue VARCHAR);\
		COPY ursa_annotations FROM STDIN DELIMITER '	';\
		INSERT INTO sample_term \
			SELECT sample.id, term.id, 1 \
			FROM ursa_annotations ursa \
			INNER JOIN sample \
			ON sample.accession=ursa.gsm \
			INNER JOIN term \
			ON term.name=ursa.tissue \
			INNER JOIN ontology \
			ON ontology.id=term.ontology_id \
			WHERE ontology.namespace='BTO';"

.PHONY: load-base load-go load-brenda

data/gene_ontology_ext.obo:
	wget -P data http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo

data/taxdump.tar.gz: 
	wget -P data/ $(NCBI_FTP)/pub/taxonomy/taxdump.tar.gz

data/gene_info.gz:
	wget -P data/ $(NCBI_FTP)/gene/DATA/gene_info.gz

data/gene2go.gz:
	wget -P data/ $(NCBI_FTP)/gene/DATA/gene2go.gz
