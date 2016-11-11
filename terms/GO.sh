wget http://purl.obolibrary.org/obo/go.obo
egrep "^id: GO|^name:" go.obo | perl -p -e 's/^id: (.*)\s*\n/GO\t\1\t/' | perl -p -e 's/name: //' | grep "GO:" | sort > GO.terms.txt

