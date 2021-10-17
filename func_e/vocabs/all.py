from .GO import getTerms as GO_getTerms
from .KEGG import getTerms as KEGG_getTerms
from .IPR import getTerms as IPR_getTerms
import pandas as pd

def getTerms(vocabs = []):
    """
    """
    terms = pd.DataFrame(columns=['Vocabulary', 'Term', 'Name'])
    if 'GO' in vocabs:
        terms = pd.concat([terms, GO_getTerms()])
    if 'IPR' in vocabs:
        terms = pd.concat([terms, IPR_getTerms()])
    if 'KEGG' in vocabs:
        terms = pd.concat([terms, KEGG_getTerms()])
    return terms
