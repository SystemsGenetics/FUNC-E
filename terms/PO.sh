wget http://purl.obolibrary.org/obo/po.obo
egrep "^id: PO|^name:" po.obo | perl -p -e 's/^id: (.*)\s*\n/PO\t\1\t/' | perl -p -e 's/name: //' | grep "PO:" | sort > PO.terms.txt
