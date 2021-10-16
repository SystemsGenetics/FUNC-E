import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as sm
from sklearn.metrics import cohen_kappa_score
import progressbar
import sys


class FUNC_E(object):

    def __init__(self):
        self.background = None
        self.terms = pd.DataFrame(columns=['Vocabulary', 'Term', 'Name'])
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
        self.enrichment = pd.DataFrame(columns=["Module", "Term", "Name", "Mod Count", "Background Count", "Fishers pVal"])
        self.clusters = None
        self.kappa = pd.DataFrame(columns=['Feature1', 'Feature2', 'Module', 'Score'])

        # The set of features that have enriched terms.
        self.efeatures = pd.DataFrame(columns=["Feature", "Module", "Term"])


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
        terms = pd.DataFrame(columns=['Vocabulary', 'Term', 'Name'])
        for tfile in files:
            new_terms = pd.read_csv(tfile, header=None, sep="\t")
            new_terms.columns = ['Vocabulary', 'Term', 'Name']
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

    def _log(self, message, level=1, end="\n"):
        if self.verbose >= level:
            print(message, end=end, flush=True)

    def run(self):
        """
        """
        if self.isReady() is False:
            self._log("Cannot perform this step as all necessary inputs are not set.")
            return

        self._log("Counting terms in the module(s) and background...")
        self.doCounts()

        self._log("Perform enrichment analysis...")
        self.doEnrichment()

        self._log("Perform multiple testing correction...")
        self.doMTC()

        self._log("Perform Kappa analysis...")
        self.doKappa()

        if self.verbose == 1:
            print("Done")

    def doModuleEnrichment(self, module):
        """
        """
        if self.isReady() is False:
            self._log("Cannot perform this step as all necessary inputs are not set.")
            return

        # Remove any previous resutls for this module that might be in the results.
        self.enrichment.drop(index=self.enrichment.index[self.enrichment['Module'] == module], inplace=True)
        self.efeatures.drop(index=self.efeatures.index[self.efeatures['Module'] == module], inplace=True)

        modCounts = self.queryCounts.loc[self.queryCounts['Module'] == module]
        modResults = pd.DataFrame(columns=["Module", "Term", "Name", "Mod Count", "Background Count", "Fishers pVal"])

        pbar = progressbar.ProgressBar(max_value=len(modCounts['Term'].unique()))

        # Second iterate through the unique modules with counts in this vocabulary.
        total_tests = 0
        for vocab in modCounts['Vocabulary'].unique():
            modVocabCounts = modCounts.loc[modCounts['Vocabulary'] == vocab]
            # Third iterate through the unique terms with counts in this module.
            for term in modVocabCounts['Term'].unique():
                total_tests = total_tests + 1
                pbar.update(total_tests)
                n11, n21, pvalue = self._performFishersTest(term, module, vocab, modCounts, modVocabCounts)

                # If the Fisher's p-value is less than the cutoff then keep it.
                if pvalue < self.ecut:
                    name = self.terms.loc[self.terms['Term'] == term]['Name'].iloc[0]
                    modResults = modResults.append({
                      "Module": module,
                      "Term": term,
                      "Name": name,
                      "Mod Count": n11,
                      "Background Count": n21,
                      "Fishers pVal": pvalue}, ignore_index=True)

        pbar.update(total_tests)
        pbar.finish()

        # Combine the module's enriched terms with the full result set.
        self.enrichment = pd.concat([self.enrichment, modResults], ignore_index=True)

        # Create the list of genes with enriched features.
        efeatures = pd.DataFrame(self.terms2features.set_index('Term')
            .join(modResults.set_index('Term'), how="inner")[['Feature', 'Module']]
            .reset_index()
            .groupby(['Feature', 'Module'])['Term']
            .apply(list)).reset_index()

        efeatures = efeatures.set_index(['Feature','Module']).join(self.query.set_index(['Feature','Module']), how='inner').reset_index()
        self.efeatures = pd.concat([self.efeatures, efeatures], ignore_index=True)


    def doEnrichment(self):
        """
        """
        if self.isReady() is False:
            self._log("Cannot perform this step as all necessary inputs are not set.")
            return

        for module in self.modules:
            self._log("Working on module: %s" % (module))
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

    def doModuleKappa(self, module):
        """
        """
        if self.isReady() is False:
            if self.verbose == 1:
                print("Cannot perform this step as all necessary inputs are not set.")
            return


        self.kappa.drop(index=self.kappa.index[self.kappa['Module'] == module], inplace=True)

        # Initialize the dataframe that will house pairwise Kappa scores.
        scores = []
        efeatures = self.efeatures[self.efeatures['Module'] == module]['Feature'].unique()

        qfterms = pd.DataFrame(self.terms2features.set_index(['Feature'])
          .join(self.query.set_index('Feature'), how='inner')
          .reset_index()
          .groupby(['Feature','Module'])['Term']
          .apply(list)).reset_index().set_index('Feature')
        qfterms = qfterms[qfterms['Module'] == module]['Term'].to_dict()
        nf = len(qfterms.keys())
        ncomps = int((nf * (nf-1)) / 2)
        self._log("Performing {:,d} Kappa comparisons for {} features...".format(ncomps, nf))

        pbar = progressbar.ProgressBar(max_value=ncomps)

        # Iterate through the list of features that have enriched terms and
        # perform pair-wise Kappa.
        total_comps = 0
        for i in range(0, len(efeatures)):
            for j in range(i+1, len(efeatures)):

                total_comps = total_comps + 1
                pbar.update(total_comps)

                fi = efeatures[i]
                fj = efeatures[j]

                fi_terms = set(qfterms[fi])
                fj_terms = set(qfterms[fj])

                u = list(fi_terms | fj_terms)

                # If we don't have enough shared terms to meet the minimum
                # overlap then skip this comparison.
                if len(u) < self.similarity_overlap:
                    continue

                li = [x in fi_terms for x in u]
                lj = [x in fj_terms for x in u]

                # Skip pairwise comparisons that have fewer overlaps that the
                # minium required.
                overlap = sum([li[x] & lj[x] for x in range(0, len(u))])

                k = 0
                if overlap < self.similarity_overlap:
                    continue
                elif overlap == len(u):
                    k = 1
                else:
                    k = cohen_kappa_score(li, lj)

                if k >= self.similarity_threshold:
                    scores.append([fi, fj, module, k])

        pbar.update(total_comps)
        pbar.finish()
        self.kappa = pd.concat([self.kappa, pd.DataFrame(scores, columns=['Feature1', 'Feature2', 'Module', 'Score'])], ignore_index=True)

    def doKappa(self):
        """
        """
        if self.isReady() is False:
            self._log("Cannot perform this step as all necessary inputs are not set.")
            return

        for module in self.modules:
            self._log("Performing Kappa test on module: %s" % (module))
            self.doModuleKappa(module)

    def _getValidSeedGroups(self, efeatures, seeds, kappa):
        """
        """
        groups = []

        # Add the seed groups that pass tests to the initial set of groups.
        for i in range(0, len(efeatures)):
            feature = efeatures[i]

            if not (feature in seeds.keys()):
                continue

            # Start with a seed with all of the features that have a passing
            # kappa score with the current feature.
            test_seed = seeds[feature]
            test_seed.add(feature)
            test_seed = list(test_seed)

            # A group must have at least a set number of genes before it
            # can be considered a seed group.
            if len(test_seed) < self.initial_group_membership:
                continue

            # keep track of the number of pairs that are good
            good_count = 0
            # keep track of total comparisions
            total_count = 0
            # Count the number of genes that have a kappa score > the threshold
            for j in range(0, len(test_seed)):
                for k in range(j+1, len(test_seed)):
                    jk_index = test_seed[j] + '-' + test_seed[k]
                    if ((jk_index in kappa.keys()) and (kappa[jk_index] > self.similarity_threshold)):
                        good_count = good_count + 1
                    total_count = total_count + 1

            # if the genes in the seed group have a set percentage (e.g. 50%) of genes
            # that have high quality kappa scores with all other genes then let's keep this
            # cluster
            if good_count / total_count >= self.percent_similarity:
                groups.append(test_seed)

        return groups

    def _mergeGroups(self, groups):
        """
        """
        print(len(groups))
        new_groups = []
        for i in range(0, len(groups)):
            for j in range(i+1, len(groups)):
                num_shared = len(set(groups[i]) & set(groups[j]))
                perci = num_shared / len(groups[i])
                percj = num_shared / len(groups[j])
                if (perci >= self.multiple_linkage_threshold) or (percj >= self.multiple_linkage_threshold):
                   new_groups.append(list(set(groups[i]) | set(groups[j])))

        if len(new_groups) > 0:
            return self._mergeGroups(new_groups)
        else:
            return new_groups


    def doModuleClustering(self, module):
        """
        """
        # Get the enriched features for this module
        efeatures = self.efeatures[self.efeatures['Module'] == module]['Feature'].sort_values().unique()

        # Generate the starting seeds by creating a dictionary of feature names
        # with the values being all of the other features with which they
        # have a meaningful kapp score.
        f1 = self.kappa.groupby('Feature1')['Feature2'].apply(list).reset_index()
        f2 = self.kappa.groupby('Feature2')['Feature1'].apply(list).reset_index()
        f1.columns = ['Feature', 'Matches']
        f2.columns = ['Feature', 'Matches']
        seeds = pd.concat([f1, f2]).groupby('Feature')['Matches'].apply(list).apply(lambda x: set(np.sort(np.unique(np.concatenate(x))))).to_dict()

        # Index the Kappa dataframe for easy lookup
        kappa =self.kappa.copy()
        kappa.index = kappa.apply(lambda x: "-".join(np.sort(np.array([x['Feature1'], x['Feature2']]))), axis=1)
        kappa = kappa['Score'].to_dict()

        groups = self._getValidSeedGroups(efeatures, seeds, kappa)
        clusters = self._mergeGroups(groups)

        return(groups)



    def doClustering(self):
        """
        """
        if self.isReady() is False:
            self._log("Cannot perform this step as all necessary inputs are not set.")
            return

        for module in self.modules:
            self._log("Performing clustering on module: %s" % (module))
            self.doModuleClustering(module)

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
        self._log("\nFisher's Test", 2)
        self._log("Module: %s" % (module), 2)
        self._log("                  Term: %s" % (term), 2)
        self._log("                    Yes      No", 2)
        self._log("               ---------------------", 2)
        self._log("In Module     | %8d | %8d | %8d" % (n11,n12,n1p), 2)
        self._log("              |----------|----------|", 2)
        self._log("In Background | %8d | %8d | %8d" % (n21,n22,n2p), 2)
        self._log("               ---------------------", 2)
        self._log("                %8d   %8d   %8d" % (np1,np2,npp), 2)
        self._log("p-value: %f" % (pvalue), 2)

        return [n11, n21, pvalue]
