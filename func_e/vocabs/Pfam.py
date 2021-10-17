import requests
import pandas as pd

def getTerms():
    return
    
    url = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz'
    r = requests.get(url, allow_redirects=True)

    terms_list = []
    in_header = True
    for line in r.content.splitlines():
        if in_header:
            in_header = False
            continue
        cols = line.decode("utf-8").split("\t")
        terms_list.append(['Pfam', cols[1], cols[4]])
    terms = pd.DataFrame(terms_list, columns=['Vocab', 'Term', 'Name'])
    return terms
