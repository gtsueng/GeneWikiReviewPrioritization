from datetime import datetime
from datetime import timedelta
import time
from Bio import Entrez
from Bio import Medline
import pandas as pd
from pandas import read_csv
import os
import re
import pickle
import lzma
import pathlib

Entrez.email = os.environ['useremail']

###############################################################################
## This module takes a list of entrez gene ids, looks up PMIDs associated with
## each gene, obtains the authors of each PMID along with their affiliations
## The returned dataframe can be processed for top authors and more
## Note that Entrez caps requests at 3/second without API key
## With API key, requests are capped at 10/second
## To be on the safe side, this module voluntarily throttles (sleeps) to
## 2 requests/second (sleep is half a second)
###############################################################################
def retrieve_detailed_pubs_by_gene(genelist):
    timestart = datetime.now().time()
    genefailures = []
    pmid_failures = []
    pmid_author_fail = []
    author_df = pd.DataFrame(columns = ["AU", "FullName","AuthorDetails","pmid","publish_date"])
    PublicationDetails = []
    print(timestart, 'obtaining publication details and authors for each gene.')
    for geneid in genelist:
        print ('fetching pmids for: '+str(geneid))
        try: 
            record = Entrez.read(Entrez.elink(dbfrom="gene", id=geneid))
            PMIDList=[] ##creates an empty list to store PMIDs   
            for link in record[0]["LinkSetDb"][0]["Link"] : ##retrieves each PMID stored in the record associated with the gene ID
                PMIDList.append(link["Id"]) ##stores PMID's linked to the gene into the list
            if len(PMIDList) > 30:
                for PMID in PMIDList: #iterates through the PMID list
                    try:
                        handle = Entrez.efetch(db="pubmed", id=int(PMID), rettype="medline", retmode="text")
                        records = Medline.parse(handle) ##parses pubmed entry for that ID and records the author
                        for record in records:
                            try:
                                PublicationDate = record.get("DP","?") #writes the publication date 
                                pubdate_type = "DP"
                            except:
                                PublicationDate = record.get("EDAT","?") #writes the initial Entrez record submission date 
                                pubdate_type = "EDAT"
                            PublicationDetails.append({'geneid':str(geneid),'pmid':PMID,'PublicationDate':PublicationDate,"PubDateType":pubdate_type})
                            try:
                                AuthorSet = record.get("AU","?") #writes the record to a list called AuthorSet                     
                                FullAuthorSet = record.get("FAU","?") #writes the record to a list called AuthorSet
                                AuthorDetails = record.get("AD","?") #writes the record to a list called AuthorDetails
                                tmp_df = pd.DataFrame({'AU':pd.Series(AuthorSet), "FullName":pd.Series(FullAuthorSet),'AuthorDetails': pd.Series(AuthorDetails)})
                                tmp_df['pmid']=PMID
                                tmp_df['publish_date']=PublicationDate
                                author_df = pd.concat((author_df,tmp_df),ignore_index=True)
                            except:
                                pmid_author_fail.append(PMID)
                        time.sleep(0.5)
                    except:
                        #print('bad pmid: ', PMID)
                        pmid_failures.append({'geneid':geneid,'pmid':PMID})
            else:
                genefailures.append(geneid)
        except:
            genefailures.append(geneid)
        
    PublicationDetailsDF = pd.DataFrame(PublicationDetails)
    PMIDfailsDF = pd.DataFrame(pmid_failures)
    timeend = datetime.now().time()
    print(timeend)
    return(PublicationDetailsDF, author_df, genefailures, pmid_author_fail, PMIDfailsDF)


################################################################################
## This uses a list of PMIDs to pull author information
################################################################################

#PMIDList = [23039619,29390967,31363486,30951672] ##Unit test

def retrieve_authors_by_pmids(PMIDList):
    print(datetime.datetime.now().time())
    author_df = pd.DataFrame(columns = ["AU", "FullName","AuthorDetails","pmid","publish_date"])
    PublicationDetails = []
    PMIDFails = []
    for PMID in PMIDList: #iterates through the PMID list
        try:
            #print('fetching authors for: '+str(PMID))
            handle = Entrez.efetch(db="pubmed", id=PMID, rettype="medline", retmode="text")
            records = Medline.parse(handle) ##parses pubmed entry for that ID and records the author
            for record in records:
                AuthorSet = record.get("AU","?") #writes the record to a list called AuthorSet
                FullAuthorSet = record.get("FAU","?") #writes the record to a list called AuthorSet
                try:
                    AuthorDetails = record.get("AD","?") #writes the record to a list called AuthorDetails
                except:
                    AuthorDetails = ['No details']
                try:
                    PublicationDate = record.get("DP","?") #writes the initial Entrez record submission date
                    PublicationDetails.append({'pmid':PMID, 'PubDateType':'PD','PubDate':PublicationDate})
                except:
                    PublicationDate = record.get("EDAT","?")
                    PublicationDetails.append({'pmid':PMID, 'PubDateType':'EDAT','PubDate':PublicationDate})
                #print(len(AuthorSet),len(AuthorDetails),AuthorSet[0],AuthorDetails[0])
                tmp_df = pd.DataFrame({'AU':pd.Series(AuthorSet), "FullName":pd.Series(FullAuthorSet),'AuthorDetails': pd.Series(AuthorDetails)})
                tmp_df['pmid']=PMID
                tmp_df['PubDate']=PublicationDate
                author_df = pd.concat((author_df,tmp_df))
        except:
            PMIDFails.append(PMID)
            print("pmid not found: ",PMID)

    PublicationDF = pd.DataFrame(PublicationDetails)
    print(datetime.datetime.now().time())
    return(PublicationDF,author_df,PMIDFails)


################################################################################
## This module parses email addresses from the AD/Affiliation field if available
## input is the author dataframe resulting from this module:
## retrieve_detailed_pubs_by_gene
## Note that Articles prior to: October 1, 2013, will have an affiliation ONLY 
## for the first author. Details on that here: 
## https://www.nlm.nih.gov/pubs/techbull/so13/brief/so13_author_affiliations.html
################################################################################

def parse_out_emails(author_df):
    author_df.reset_index(inplace=True)
    try:
        author_df.drop(['level_0'],axis=1,inplace=True)
    except Exception: 
        pass
    try:
        author_df.drop(['index'],axis=1,inplace=True)
    except Exception: 
        pass
    author_df['email'] = author_df['AuthorDetails'].str.extract(r'([^@|\s]+@[^@]+\.[^@|\s]+)')
    author_emails = author_df.loc[author_df['email'].notnull()]
    author_emails_dots = author_emails.loc[author_emails['email'].str[-1]=="."]
    author_emails_dots['email'] = author_emails_dots['email'].str[:-1]
    author_emails.update(author_emails_dots)
    author_df.update(author_emails)
    return(author_df)


###############################################################################
## This module takes the resulting dataframes from the previous module
## merges them and then gets a count of the number of publications an author
## has for each gene. Note that the optional parameter, method allows for the
## counts to be based on the author "AU" or full name of the author "FullName"
## The default is the FullName
###############################################################################
def get_top_authors_from_dfs(PublicationDetailsDF, author_df, method="FullName"):
    all_merged_df = author_df_deets.merge(PublicationDetailsDF, on='pmid',how='left')
    all_merged_df.drop_duplicates(keep='first',inplace=True)
    top_authors_per_gene = all_merged_df.groupby(['geneid',method,'FullName']).size().reset_index(name='pubcounts')
    top_authors_per_gene.sort_values(by=['geneid','pubcounts'],ascending=[True,False],inplace=True)
    return(top_authors_per_gene)


## Pull up all pmids per genes, get author details, pull out available email addresses
## Perform groupby counts to get top contributing authors per gene
        
def get_authors(genelist,datapath,test=False):
    if test==True:
        genelist = [439921,55768] ## for unit test     
    PublicationDetailsDF, author_df, genefailures, pmid_author_fail, PMIDfailsDF = retrieve_detailed_pubs_by_gene(genelist)
    author_df_deets = parse_out_emails(author_df)
    PublicationDetailsDF.to_csv(os.path.join(datapath,'PublicationDetailsDF.tsv'),sep='\t',header=True)
    #author_df.to_csv(os.path.join(datapath,'author_df.tsv'),sep='\t',header=True)
    author_pickle = pickle.dumps(author_df)
    with lzma.open(os.path.join(datapath,"author_df.xz"), "w") as f:
        f.write(author_pickle)
    #author_df_deets.to_csv(os.path.join(datapath,'author_df_deets.tsv'),sep='\t',header=True)
    PMIDfailsDF.to_csv(os.path.join(datapath,'PMIDfailsDF.tsv'),sep='\t',header=True)
    with open(datapath+'genefailures.txt','w') as outwrite:
        for eachgene in genefailures:
            outwrite.write(str(eachgene)+'\n')
        outwrite.close()

def deal_with_failures(datapath,hasfailures = True):
    oldPublicationDetailsDF = read_csv(os.path.join(datapath,'PublicationDetailsDF.tsv'),delimiter='\t',header=0,index_col=0)
    oldauthor_df = read_csv(os.path.join(datapath,'author_df.tsv'),delimiter='\t',header=0,index_col=0)
    oldauthor_df_deets = read_csv(os.path.join(datapath,'author_df_deets.tsv'),delimiter='\t',header=0,index_col=0)
    oldPMIDfailsDF = read_csv(os.path.join(datapath,'PMIDfailsDF.tsv'),delimiter='\t',header=0,index_col=0)
    i=0
    while hasfailures == True:
        with open(datapath+'genefailures.txt','r') as infile:
            tmpgenelist = []
            for line in infile:
                tmpgenelist.append(line.strip())
            infile.close()
        tmpPublicationDetailsDF, tmpauthor_df, tmpgenefailures, tmppmid_author_fail, tmpPMIDfailsDF = retrieve_detailed_pubs_by_gene(tmpgenelist)
        tmpauthor_df_deets = parse_out_emails(author_df)
        PublicationDetailsDF = pd.concat((oldPublicationDetailsDF,tmpPublicationDetailsDF),ignore_index=True)
        PublicationDetailsDF.to_csv(os.path.join(datapath,'PublicationDetailsDF.tsv'),sep='\t',header=True)
        author_df = pd.concat((oldauthor_df,tmpauthor_df),ignore_index=True)
        author_df.to_csv(os.path.join(datapath,'author_df.tsv'),sep='\t',header=True)
        author_df_deets = pd.concat((oldauthor_df_deets,tmpauthor_df_deets),ignore_index=True)
        author_df_deets.to_csv(os.path.join(datapath,'author_df_deets.tsv'),sep='\t',header=True)
        PMIDfailsDF = pd.concat((oldPMIDfailsDF,tmpPMIDfailsDF),ignore_index=True)
        PMIDfailsDF.to_csv(os.path.join(datapath,'PMIDfailsDF.tsv'),sep='\t',header=True)
        with open(datapath+'genefailures.txt','w') as outwrite:
            for eachgene in tmpgenefailures:
                outwrite.write(str(eachgene)+'\n')
            outwrite.close()
        if len(tmpgenefailures) < 2:
            hasfailures = False
        i=i+1
    return(i)


## Pull top authors for high priority genes
script_path = pathlib.Path(__file__).parent.absolute()
datapath = os.path.join(script_path,'data/')
resultpath = os.path.join(script_path,'results/')
genefile = 'priority_by_size.tsv'

prioritylist = read_csv(os.path.join(resultpath,genefile),delimiter='\t',header=0,index_col=0)
genelist = prioritylist['geneID'].unique().tolist()

get_authors(genelist,datapath,test=False)




