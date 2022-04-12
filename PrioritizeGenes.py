import os
import pandas as pd
from pandas import read_csv
import requests
import pathlib



def load_authordf(datapath):
    import lzma
    import pickle
    with lzma.open(os.path.join(datapath,"author_df.xz")) as f:
        authorpickle = f.read()
    authordf = pickle.loads(authorpickle)
    return(authordf)


def get_clean_authors(datapath):
    author_df = load_authordf(datapath)
    author_sum = author_df.groupby('AU').size().reset_index(name='count')
    no_single_authors = author_sum.loc[author_sum['count']>1].copy()
    less_details = author_df.drop(['AuthorDetails','publish_date'],axis=1).copy()
    clean_authors = less_details.loc[less_details['AU'].isin(no_single_authors['AU'].tolist())]
    return(clean_authors)


def get_gene_pmid_table(datapath):
    pub_details = read_csv(os.path.join(datapath,'PublicationDetailsDF.tsv'), delimiter='\t',index_col=0,header=0)
    gene_pmid = pub_details[['geneid','pmid']].copy()
    gene_pmid.drop_duplicates(keep='first',inplace=True)
    return(gene_pmid)


def generate_pub_summary(datapath):
    pub_details = read_csv(os.path.join(datapath,'PublicationDetailsDF.tsv'), delimiter='\t',index_col=0,header=0)
    pub_details['year'] = pub_details['PublicationDate'].str.extract(r'(\d\d\d\d)')
    pub_details.drop('PubDateType',axis=1,inplace=True)
    pub_details.drop_duplicates(keep='first',inplace=True)
    pub_details['year'] = pub_details['year'].fillna(0).astype(int)
    pub_frequency = pub_details.groupby('geneid').size().reset_index(name='pubcount')
    median_year = pub_details.groupby('geneid')['year'].median().reset_index(name='median_pub_year')
    median_year['median_pub_year'] = median_year['median_pub_year'].astype(int)
    max_year = pub_details.groupby('geneid')['year'].max().reset_index(name='max_pub_year')
    max_year['max_pub_year'] = max_year['max_pub_year'].astype(int)
    pub_sum = pub_frequency.merge(median_year.merge(max_year,on='geneid',how='left'),on='geneid',how='left').copy()
    return(pub_sum)


def merge_and_filter_results(datapath,resultpath,min_pubcount=30,min_pagelength=200):
    priority_by_size = read_csv(os.path.join(resultpath,'priority_by_size.tsv'), delimiter='\t',index_col=0,header=0)
    priority_by_size.rename(columns={'geneID':'geneid'},inplace=True)
    pub_sum = generate_pub_summary(datapath)
    gene_summary = priority_by_size.merge(pub_sum,on='geneid',how='inner')
    filtered_gene_summary = gene_summary.loc[((gene_summary['pubcount']>min_pubcount) & (gene_summary['page_length']>min_pagelength))].copy()
    filtered_gene_summary.sort_values('page_length',ascending=True,inplace=True)
    filtered_gene_summary.head(n=500).to_csv(os.path.join(resultpath,'genes_by_wiki_length.tsv'),sep='\t',header=True)
    scored_gene_summary = filtered_gene_summary.copy()
    scored_gene_summary['priority_score']= 10000/scored_gene_summary['page_length']+(scored_gene_summary['pubcount']/2800)
    scored_gene_summary.sort_values('priority_score',ascending=False,inplace=True)
    scored_gene_summary.head(n=500).to_csv(os.path.join(resultpath,'genes_by_score.tsv'),sep='\t',header=True)
    scored_gene_summary.sort_values('pubcount',ascending=False,inplace=True)
    scored_gene_summary.head(n=500).to_csv(os.path.join(resultpath,'genes_by_pubcount.tsv'),sep='\t',header=True)
    
    
def generate_author_table(datapath,resultpath):
    author_by_gene = pd.DataFrame(columns=['geneid','AU','counts','FullName','email'])
    gene_pmid = get_gene_pmid_table(datapath)
    clean_authors = get_clean_authors(datapath)
    for eachgene in gene_pmid['geneid'].unique().tolist():
        pmids = gene_pmid['pmid'].loc[gene_pmid['geneid']==eachgene]
        tmpauths = clean_authors.loc[clean_authors['pmid'].isin(pmids)]
        cleanauths = tmpauths.drop_duplicates(subset=['AU','FullName','pmid'],keep='first')
        authorlist = cleanauths.groupby(['AU']).size().reset_index(name='counts')
        to_keep = authorlist.loc[authorlist['counts']>2].copy()
        if len(to_keep)>0:
            auth_deets = to_keep.merge(tmpauths,on='AU',how='inner')
            auth_deets.drop('pmid',axis=1,inplace=True)
            auth_deets.drop_duplicates(keep='first',inplace=True)
            auth_deets['geneid']=eachgene
            auth_deets.sort_values('counts',ascending=False,inplace=True)
            author_by_gene = pd.concat((author_by_gene,auth_deets),ignore_index=True)
    author_by_gene.to_csv(os.path.join(resultpath,'potential_authors.tsv'),sep='\t',header=True)
    

#### Main Function
script_path = pathlib.Path(__file__).parent.absolute()
datapath = os.path.join(script_path,'data/')
resultpath = os.path.join(script_path,'results/')
generate_author_table(datapath,resultpath)
