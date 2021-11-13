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

            if len(cols.keys()) == 5:
                terms_list.append([cols['ID_Space'], cols['Vocabulary'], cols['Term'], cols['Name'], cols['Definition']])
            cols = {}
            m = re.search(r'^id: (GO:\d+)', line)
            cols['Term'] = m.group(1)
        if re.search(r'^name: ', line):
            m = re.search(r'^name: (.+)', line)
            cols['Name'] = m.group(1)
        if re.search(r'^namespace: ', line):
            m = re.search(r'^namespace: (.+)', line)
            cols['ID_Space'] = 'GO'
            cols['Vocabulary'] = m.group(1)
        if  re.search(r'^def: ', line):
            m = re.search(r'^def: "(.+)"', line)
            cols['Definition'] = m.group(1)
    terms = pd.DataFrame(terms_list, columns=['ID_Space', 'Vocabulary', 'Term', 'Name', 'Definition'])
    return terms
