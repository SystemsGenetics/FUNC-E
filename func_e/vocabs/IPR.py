import requests
import pandas as pd

def getTerms():
    url = 'http://ftp.ebi.ac.uk/pub/databases/interpro/entry.list'
    r = requests.get(url, allow_redirects=True)

    terms_list = []
    in_header = True
    for line in r.content.splitlines():
        if in_header:
            in_header = False
            continue
        cols = line.decode("utf-8").split("\t")
        terms_list.append(['IPR', cols[0], cols[2]])
    terms = pd.DataFrame(terms_list, columns=['Vocabulary', 'Term', 'Name'])
    return terms
