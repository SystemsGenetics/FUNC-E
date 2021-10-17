import requests
import pandas as pd

def getOrthologs():
    """
    """
    url = 'http://rest.kegg.jp/list/ko'
    r = requests.get(url, allow_redirects=True)

    terms_list = []
    for line in r.content.splitlines():
        cols = line.decode("utf-8").split("\t")
        terms_list.append(['KEGG', cols[0], cols[1]])
    terms = pd.DataFrame(terms_list, columns=['Vocab', 'Term', 'Name'])
    terms['Term'] = terms['Term'].str.replace(r'ko:','', regex=True)

    return terms

def getPathways():
    """
    """
    url = 'http://rest.kegg.jp/list/pathway'
    r = requests.get(url, allow_redirects=True)

    terms_list = []
    for line in r.content.splitlines():
        cols = line.decode("utf-8").split("\t")
        terms_list.append(['KEGG', cols[0], cols[1]])
    terms = pd.DataFrame(terms_list, columns=['Vocab', 'Term', 'Name'])
    terms['Term'] = terms['Term'].str.replace(r'path:map','ko', regex=True)
    return terms

def getModules():
    """
    """
    url = 'http://rest.kegg.jp/list/md'
    r = requests.get(url, allow_redirects=True)

    terms_list = []
    for line in r.content.splitlines():
        cols = line.decode("utf-8").split("\t")
        terms_list.append(['KEGG', cols[0], cols[1]])
    terms = pd.DataFrame(terms_list, columns=['Vocab', 'Term', 'Name'])
    terms['Term'] = terms['Term'].str.replace(r'md:','', regex=True)
    return terms

def getTerms():
    modules = getModules()
    orthologs = getOrthologs()
    pathways = getPathways()

    return pd.concat([pathways, orthologs, modules])
