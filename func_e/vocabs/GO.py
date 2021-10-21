import requests
import re
import pandas as pd

def getTerms():
    url = 'http://purl.obolibrary.org/obo/go.obo'
    r = requests.get(url, allow_redirects=True)

    terms_list = []
    cols = {}
    for line in r.content.splitlines():
        line = line.decode("utf-8")
        if re.search(r'^id: GO', line):
            if len(cols.keys()) == 3:
                terms_list.append([cols['Vocabulary'], cols['Term'], cols['Name']])
            cols = {}
            m = re.search(r'^id: (GO:\d+)', line)
            cols['Term'] = m.group(1)
        if re.search(r'^name: ', line):
            m = re.search(r'^name: (.+)', line)
            cols['Name'] = m.group(1)
        if re.search(r'^namespace: ', line):
            m = re.search(r'^namespace: (.+)', line)
            cols['Vocabulary'] = 'GO' #m.group(1)
    terms = pd.DataFrame(terms_list, columns=['Vocabulary', 'Term', 'Name'])
    return terms
