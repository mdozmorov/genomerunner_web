#!/usr/bin/env bash

###########################
# User-modifiable variables
###########################

DB=pubmed

#############
# Main script
#############

usage() {
    cat <<EOF
USAGE: $0 path/to/medline/gz/dir
EOF
}

[ $# -ne 1 ] && {
    usage
    exit 1
}

echo "Creating database and (re)loading schema ..." 1>&2

hgsql -e "CREATE DATABASE IF NOT EXISTS $DB;"
hgsql $DB < schema.sql

for file in $1/*.xml.gz; do 
    echo Loading $file ... 1>&2
    # FIXME (?) Use REPLACE or the default (IGNORE) ?
    zcat $file | xsltproc --novalid article.xsl - \
        | grep -v "^$" | hgsql $DB -e \
            "LOAD DATA LOCAL INFILE '/dev/stdin' REPLACE INTO TABLE article"
done
