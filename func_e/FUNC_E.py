import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as sm


class FUNC_E(object):

    def __init__(self):
        self.background = None
        self.terms = pd.DataFrame(columns=['Vocabulary', 'Term', 'Definition'])
        self.modules = None
        self.query = None
        self.terms2features = pd.DataFrame(columns=['Feature', 'Term'])

        self.bgCounts = None
        self.queryCounts = None

        # If set to a non-zero value then logging will occur
        self.verbose = 0

        # Enrichment setting defaults
        self.ecut = 0.05

        # Clustering setting defaults
        self.similarity_overlap = 3
        self.percent_similarity = 0.50
        self.initial_group_membership = 3
        self.multiple_linkage_threshold = 0.50
        self.final_group_membership = 3
        self.similarity_threshold = 0.50

        # Dataframes containing results
        self.enrichment = pd.DataFrame(columns=["Module", "Term", "Definition", "Mod Count", "Background Count", "Fishers pVal"])
        self.clusters = None

    def setVerbosity(self, level=0):
        self.verbose = level

    def setEnrichmentSettings(self, settings):
        if settings['ecut']:
            self.ecut = settings['ecut']

    def setClusteringSettings(self, settings):
        if settings['similarity_overlap']:
            self.similarity_overlap = settings['similarity_overlap']
        if settings['percent_similarity']:
            self.percent_similarity = settings['percent_similarity']
        if settings['initial_group_membership']:
            self.initial_group_membership = settings['initial_group_membership']
        if settings['multiple_linkage_threshold']:
            self.multiple_linkage_threshold = settings['multiple_linkage_threshold']
        if settings['final_group_membership']:
            self.final_group_membership = settings['final_group_membership']

    def applyClusteringPreset(self, preset):
        """
        """
        if preset == "lowest":
            self.similarity_threshold = 0.20
        if preset == "low":
            self.similarity_threshold = 0.35
        if preset == "medium":
            self.similarity_threshold = 0.50
        if preset == "high":
            self.similarity_threshold = 0.85
        if preset == "highest":
            self.similarity_threshold = 1.00

    def setBackground(self, background):
        # TODO: add some quality checks to make sure the
        # background dataframe is appropriate.
        self.background = background

    def importBackgroundFile(self, file):
        background = pd.read_csv(file, header=None)
        background.columns = ['Feature']
        background.set_index('Feature')
        self.setBackground(background)

    def setQuery(self, query):
        self.query = query

    def importQueryFile(self, file):
        query = pd.read_csv(file, header=None, sep="\t")
        if len(query.columns) == 1:
            query.columns = ['Feature']
            query['Module'] = 'module0'
        else:
            query.columns = ['Feature', 'Module']
        self.setQuery(query)

    def setTerms(self, terms):
        self.terms = terms

    def importTermsFiles(self, files):
        terms = pd.DataFrame(columns=['Vocabulary', 'Term', 'Definition'])
        for tfile in files:
            new_terms = pd.read_csv(tfile, header=None, sep="\t")
            new_terms.columns = ['Vocabulary', 'Term', 'Definition']
            terms = pd.concat([terms, new_terms])
        self.setTerms(terms)

    def setTerms2Features(self, terms2features):
        self.terms2features = terms2features

    def importTerms2FeaturesFiles(self, files):
        """
        """
        terms2features = pd.DataFrame(columns=['Feature', 'Term'])
        for t2ffile in files:
            new_terms2f = pd.read_csv(t2ffile, header=None, sep="\t")
            new_terms2f.columns = ['Feature', 'Term']
            terms2features = pd.concat([terms2features, new_terms2f])
        terms2features = terms2features.set_index('Feature', drop=False)
        self.setTerms2Features(terms2features)

    def doCounts(self):
        """
        Counts the number of times terms are present in the background and modules.

        :param background:  the Pandas dataframe containing the background
        features.

        :param query: the Pandas dataframe containing the query features and the
        modules they belong to.

        :param terms: the Pandas dataframe containing the list of vocabulary terms.

        :param terms2features: the Pandas dataframe mapping the vocabulary terms
        to the background features.

        :return:  An array of two Pandas dataframes, one containing the
        background counts and the other containing the query module counts.
        """

        # First merge the terms and term to feature mapping.
        t2f_full = self.terms2features.set_index('Term').join(self.terms.set_index('Term'), on='Term', how="left")
        t2f_full = t2f_full.reset_index()

        # Count the background terms.
        bg2terms = t2f_full.set_index('Feature').join(self.background.set_index('Feature'), on='Feature', how="left")
        bg2terms = bg2terms.reset_index()
        self.bgCounts = bg2terms.groupby(['Vocabulary', 'Term']).nunique()
        self.bgCounts = self.bgCounts['Feature'].reset_index()

        # Count the terms in the query per module.
        queryTerms = t2f_full.set_index('Feature').join(self.query.set_index('Feature'), on='Feature', how="left")
        queryTerms = queryTerms.reset_index()
        self.queryCounts = queryTerms.groupby(['Module','Vocabulary','Term']).nunique()
        self.queryCounts = self.queryCounts['Feature'].reset_index()
        self.modules = self.queryCounts['Module'].unique()

    def isReady(self):
        """
        """
        if self.query is None:
            return False
        if self.background is None:
            return False
        if self.terms.shape[0] == 0:
            return False
        if self.terms2features.shape[0] == 0:
            return False
        return True

    def doModuleEnrichment(self, module):
        """
        """
        if self.isReady() is False:
            if self.verbose == 1:
                print("Cannot perform this step as all necessary inputs are not set.")
            return

        modCounts = self.queryCounts.loc[self.queryCounts['Module'] == module]
        modResults = pd.DataFrame(columns=["Module", "Term", "Definition", "Mod Count", "Background Count", "Fishers pVal"])

        # Second iterate through the unique modules with counts in this vocabulary.
        for vocab in modCounts['Vocabulary'].unique():
            modVocabCounts = modCounts.loc[modCounts['Vocabulary'] == vocab]
            # Third iterate through the unique terms with counts in this module.
            for term in modVocabCounts['Term'].unique():
                n11, n21, pvalue = self._performFishersTest(term, module, vocab, modCounts, modVocabCounts)

                # If the Fisher's p-value is less than the cutoff then keep it.
                if pvalue < self.ecut:
                    definition = self.terms.loc[self.terms['Term'] == term]['Definition'].iloc[0]
                    modResults = modResults.append({"Module": module, "Term": term, "Definition": definition, "Mod Count": n11, "Background Count": n21, "Fishers pVal": pvalue}, ignore_index=True)
        self.enrichment = pd.concat([self.enrichment, modResults], ignore_index=True)

    def doEnrichment(self):
        """
        """
        if self.isReady() is False:
            if self.verbose == 1:
                print("Cannot perform this step as all necessary inputs are not set.")
            return

        for module in self.modules:
            if self.verbose > 0:
                print("Performing Fisher's Tests on module: %s" % (module))
            self.doModuleEnrichment(module)

    def doMTC(self):
        """
        Apply multiple testing correction using Bonferroni and Benjamini-Hochberg
        on a per-module basis.
        """
        if self.isReady() is False:
            if self.verbose == 1:
                print("Cannot perform this step as all necessary inputs are not set.")
            return

        bonferroni = [None, None]       # Default length-two list for scope
        benjamini = [None, None]
        if len(self.enrichment["Fishers pVal"]) > 0:     # some terms are significant by ecut standard
            bonferroni = sm.multipletests(self.enrichment["Fishers pVal"], method='bonferroni')
            benjamini = sm.multipletests(self.enrichment["Fishers pVal"], method='fdr_bh')
        else:
            bonferroni = ["Not enough significant terms", "Not enough significant terms"]   # message in result if insufficient terms
            benjamini = ["Not enough significant terms", "Not enough significant terms"]
        self.enrichment['Bonferroni'] = bonferroni[1]
        self.enrichment['Benjamini'] = benjamini[1]

    def doModuleClustering(self, module):
        if self.isReady() is False:
            if self.verbose == 1:
                print("Cannot perform this step as all necessary inputs are not set.")
            return

        # Initialize the dataframe that will house pairwise Kappa scores.
        kappaResults = pd.DataFrame(columns=['Feature1', 'Feature2', 'Score'])

        modResults = self.enrichment[self.enrichment['Module'] == module]

        # Get the list of features that have enriched terms.
        qModule = self.query.loc[self.query['Module'] == module]
        qModTerms = qModule.set_index('Feature').join(self.terms2features.set_index('Feature'), on="Feature", how="left", lsuffix="_q", rsuffix="_t2f")
        qModTerms = qModTerms.reset_index()
        efeatures = modResults.set_index('Term').join(qModTerms.set_index('Term'), on="Term", lsuffix="_res", rsuffix="_q")
        efeatures = efeatures['Feature'].unique()
        efeatures.sort()

        # Iterate through the list of features that have enriched terms and
        # perform pair-wise Kappa.
        for i in range(0, len(efeatures)):
            if self.verbose > 0:
                print("  Working on feature %d of %d" % (i, len(efeatures)))
            for j in range(i+1, len(efeatures)):
                pass
                k = self._calculateKappa(efeatures[i], efeatures[j])
                if k >= self.similarity_threshold:
                    kappaResults = kappaResults.append({'Feature1': efeatures[i], 'Feature2': efeatures[j], 'Score': k}, ignore_index=True)
                    if self.verbose > 2:
                        print("%d of %d, %s vs %s: %f" % (j, len(efeatures), efeatures[i], efeatures[j], k))

        kappaResults.index = pd.MultiIndex.from_frame(kappaResults[['Feature1', 'Feature2']])

    def doClustering(self):
        """
        """
        if self.isReady() is False:
            if self.verbose == 1:
                print("Cannot perform this step as all necessary inputs are not set.")
            return

        for module in self.modules:
            if self.verbose > 0:
                print("Performing Clustering on module: %s" % (module))
            self.doClustering(module)

    def _performFishersTest(self, term, module, vocab, modCounts, modVocabCounts):
        """
        Performs the Fisher's Exact Test to see if a term is enriched in a module.

        :param term: the name of the term to test

        :param module:  the name of the module where counts should be retrieved.

        :param vocabulary:  the vocabulary name to which the term belongs.

        :param modCounts: a Pandas Dataframe containing the count of all terms
        assigned to features in a module.

        :param modVocabCounts: a Pandas dataframe containing the count of all terms
        from the vocabulary assigned to features in the module.

        :param bgCounts:  a Pandas dataframe containing the count of all terms
        assigned to features in the background.

        :return:  A list containing the number of terms in the module, the
        number of terms in the background and the Fisher's p-value.
        """
        #  Contigency matrix for each term in a module:
        #
        #                     Yes       No     Totals
        #                  ------------------
        #  In Module       |  n11   |   n12  | n1p
        #  In Background   |  n21   |   n22  | n2p
        #                  -----------------
        #  Totals           np1      np2      npp
        #
        n11 = modCounts.loc[modCounts['Term'] == term]['Feature'].iloc[0]
        n21 = self.bgCounts.loc[self.bgCounts['Term'] == term]['Feature'].iloc[0]
        n1p = modVocabCounts['Feature'].sum()
        n2p = self.bgCounts.loc[self.bgCounts["Vocabulary"] == vocab]['Feature'].sum()
        n12 = n1p - n11;
        n22 = n2p - n21;
        np1 = n11 + n21;
        np2 = n12 + n22;
        npp = np1 + np2;
        oddsratio, pvalue = stats.fisher_exact([[n11, n12], [n21, n22]], alternative="greater")
        if self.verbose > 1:
            print("\nFisher's Test")
            print("Module: %s" % (module))
            print("                  Term: %s" % (term))
            print("                    Yes      No")
            print("               ---------------------")
            print("In Module     | %8d | %8d | %8d" % (n11,n12,n1p))
            print("              |----------|----------|")
            print("In Background | %8d | %8d | %8d" % (n21,n22,n2p))
            print("               ---------------------")
            print("                %8d   %8d   %8d" % (np1,np2,npp))
            print("p-value: %f" % (pvalue))

        return [n11, n21, pvalue]

    def _calculateKappa(self, feature1, feature2, similarity_overlap, similarity_threshold):
        """

        """
        # Get the lsit of terms assigned to each feature and join the lists.
        # This joining will allow us to see which terms are in common.
        i = set(self.terms2features.loc[feature1, 'Term'])
        j = set(self.terms2features.loc[feature2, 'Term'])

        # The contigency matrix is used for calculating Kappa.
        #
        #                 Gene i
        #
        #  G       |   In  |  Not  |  total
        #  e     --|-------|-------|-------
        #  n    In |  c11  |  c10  |  c1_
        #  e     --|-------|-------|-------
        #      Not |  c01  |  c00  |  c0_
        #  j     --|-------|-------|-------
        #    total |  c_1  |  c_0  |  tab
        #
        #  c11 = number of terms in common between gene i and gene j
        #  c10 = number of terms in gene j but not in gene i
        #  c01 = number of terms in gene i but not in gene j
        #  c00 = number of terms in neither gene i nor gene j
        #
        c11 = len(i.intersection(j))
        c10 = len(i.difference(j))
        c01 = len(j.difference(i))
        c00 = self.bgCounts['Feature'].sum() - (c01 + c10 + c11)
        c_1 = c11 + c01
        c_0 = c10 + c00
        c0_ = c01 + c00
        c1_ = c11 + c10
        tab = c1_ + c0_

        # Don't perform kappa stats on genes that share less than
        # similarity_overlap number of terms.
        if c11 < similarity_overlap:
            return 0

        # Make sure our counts are in agreement
        if c1_ + c0_ != c_1 + c_0:
           print("Kappa cannot be calculated because the number of agreements in gene i and gene j are not equal")
           exit(1)

        # Calculate the kappa score
        oa = (c11 + c00) / tab
        ca = (c_1 * c1_ + c_0 * c0_) / (tab * tab)

        # Skip this if the chance agreement == 1
        if(ca == 1):
            return 0

        # Calculate the Kappa score.
        k = (oa - ca) / (1 - ca)

        # If verbosify level is sufficient then let's print out the
        # successful Kappa tests.
        if self.verbose > 1:
            if k >= similarity_threshold:
                print("\nKappa Stats")
                print("  Comparison: %s vs %s" % (feature1, feature2))
                print("            In i      Not i")
                print("         ---------------------")
                print(" In J   | %8d | %8d | %8d" % (c11, c01, c_1))
                print("        |----------|----------|")
                print(" Not J  | %8d | %8d | %8d" % (c10, c00, c_0))
                print("         ---------------------")
                print("          %8d   %8d   %8d" % (c1_, c0_, tab))
                print("   Observed:  %f" % (oa))
                print("   Chance:    %f" % (ca))
                print("   Kappa:     %f" % (k))
                if self.verbose > 2:
                    print(i)
                    print(j)

        return k
