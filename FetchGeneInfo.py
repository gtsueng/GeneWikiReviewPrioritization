#### Import libraries
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import pandas as pd
import os
import json
import urllib.parse
import urllib.request
import time
import mwclient as mw
from datetime import datetime
import pathlib

###############################################################################
## Request nicely
###############################################################################

DEFAULT_TIMEOUT = 5 # seconds

class TimeoutHTTPAdapter(HTTPAdapter):
    def __init__(self, *args, **kwargs):
        self.timeout = DEFAULT_TIMEOUT
        if "timeout" in kwargs:
            self.timeout = kwargs["timeout"]
            del kwargs["timeout"]
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        timeout = kwargs.get("timeout")
        if timeout is None:
            kwargs["timeout"] = self.timeout
        return super().send(request, **kwargs)

## Set time outs, backoff, retries
httprequests = requests.Session()
retry_strategy = Retry(
    total=3,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
    method_whitelist=["HEAD", "GET", "OPTIONS"] ## Note this method is deprecated and replaced with `allowed_methods` for newer releases of requests library
    #allowed_methods=["HEAD", "GET", "OPTIONS"] ## Note this method is deprecated and replaced with `allowed_methods` for newer releases of requests library
)
adapter = TimeoutHTTPAdapter(timeout=25,max_retries=retry_strategy)
httprequests.mount("https://", adapter)
httprequests.mount("http://", adapter)

###############################################################################
## This module uses mwclient to pull page size and edit stats on wikipedia pages  
## for each gene given a list of gene wikipedia titles
###############################################################################
def get_wiki_volume_info (mwsite,titlelist):
    print('obtaining wikipedia volume information')
    titlelist
    pageinfo=[]
    pagefails = []
    for eachpage in titlelist:
        tempdict={} #title, length/size, last_revised, last_revision_id
        try:
            checkitem = mwsite.api('query', prop='info', titles=eachpage)
            results1 = checkitem['query']['pages']
            for item in results1:
                base = str(item)
                results2 = results1[base]
                tempdict['title']=str(results2['title'])
                tempdict['page_length']=int(results2['length'])
                tempdict['last_touched']=str(results2['touched'])
                tempdict['lastrevid']=str(results2['lastrevid'])
                pageinfo.append(tempdict)               
        except:
            pagefails.append(eachpage)
            pass 
        time.sleep(1)
    return(pageinfo,pagefails)

###############################################################################
## This module uses pulls pageview data from the Media Wiki PageViews API
## More on the API here: https://wikimedia.org/api/rest_v1/#/Pageviews%20data/
## The module pulls in a parameter dictionary, and the list of wiki titles
## Parameters include:
## project: en.wikipedia.org, other wikimedia projects
## access: all-access, desktop, mobile-app, mobile-web
## agent: all-agents, user, spider, bot
## granularity: daily, monthly
###############################################################################
def get_monthly_pvs(page_view_parameters, useragent, no_missing):
    no_missing['titlelist'] = [x.replace(" ","_").replace("https://","http://").replace("http://en.wikipedia.org/wiki/","") for x in no_missing['Gene Wiki Page']]
    pginfo = []
    pgfails = []
    print('obtaining wikipedia pageview information')
    pv_api_url = "https://wikimedia.org/api/rest_v1/metrics/pageviews/per-article/en.wikipedia/"
    for eachtitle in no_missing['titlelist']:
        try:
            url = pv_api_url+pv_params['access']+pv_params['agent']+eachtitle+"/"+pv_params['granularity']+pv_params['start']+"/"+pv_params['end']
            r = httprequests.get(url, headers=useragent)
            items = r.json()
            try:
                for item in items["items"]:
                    tmpdict = {'title':item["article"], 'views':int(item["views"]), 'granularity':item['granularity'],
                               'timestamp':item["timestamp"],'access':item['access'],'agent':item['agent']}
                    pginfo.append(tmpdict)
            except:
                tmpdict = {'title':title, 'views':-1, 'granularity':"no data",
                               'timestamp':"00000000",'access':"not data",'agent':"no data"}
                pginfo.append(tmpdict)            
        except:
            pgfails.append(eachtitle)
        time.sleep(1)

    pginfodf = pd.DataFrame(pginfo)
    
    return(pginfodf, pgfails)    

#### Generate table of all human genes that have do not have a Wikipedia article
def get_genes_no_wiki(datapath):
    url = 'https://query.wikidata.org/sparql'
    query = """
        SELECT ?item ?itemLabel ?geneID ?proteinwdid
        WHERE
        {
          ?item wdt:P31 wd:Q7187 .
          ?item wdt:P703 wd:Q15978631 .
          ?item wdt:P351 ?geneID .
          ?item wdt:P688 ?proteinwdid .
          ?sitelink schema:about ?item .
          FILTER NOT EXISTS {
            ?article schema:about ?item .
            ?article schema:isPartOf <https://en.wikipedia.org/> .       }

        SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en" }
        }
    """
    r = httprequests.get(url, params = {'format': 'json', 'query': query})
    data = r.json()
    datadf = pd.DataFrame(data['results']['bindings'])
    datadf['uri'] = [x['value'] for x in datadf['item']]
    datadf['label'] = [x['value'] for x in datadf['itemLabel']]
    datadf['geneID'] = [x['value'] for x in datadf['geneID']]
    datadf['QID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['uri']]
    datadf['proteinuri'] = [x['value'] for x in datadf['proteinwdid']]
    datadf['proteinID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['proteinuri']]
    cleandata = datadf[['QID','label','geneID','proteinID']].copy()
    cleandata.drop_duplicates(keep='first', inplace=True)
    cleandata.to_csv(os.path.join(datapath,'genes_no_wiki.tsv'),sep = '\t', header=True)
    print(len(cleandata))
    
#### Generate table of all human genes, where the protein encoded by those genes do not have a Wikipedia article
def get_gene_proteins_no_wiki(datapath):
    url = 'https://query.wikidata.org/sparql'
    query = """
        SELECT ?item ?itemLabel ?geneID ?proteinwdid
        WHERE
        {
          ?item wdt:P31 wd:Q7187 .
          ?item wdt:P703 wd:Q15978631 .
          ?item wdt:P351 ?geneID .
          ?item wdt:P688 ?proteinwdid .
          ?sitelink schema:about ?item .
          FILTER NOT EXISTS {
            ?article schema:about ?proteinwdid .
            ?article schema:isPartOf <https://en.wikipedia.org/> .       }

        SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en" }
        }
    """
    r = httprequests.get(url, params = {'format': 'json', 'query': query})
    data = r.json()
    datadf = pd.DataFrame(data['results']['bindings'])
    datadf['uri'] = [x['value'] for x in datadf['item']]
    datadf['label'] = [x['value'] for x in datadf['itemLabel']]
    datadf['geneID'] = [x['value'] for x in datadf['geneID']]
    datadf['QID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['uri']]
    datadf['proteinuri'] = [x['value'] for x in datadf['proteinwdid']]
    datadf['proteinID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['proteinuri']]
    cleandata = datadf[['QID','label','geneID','proteinID']].copy()
    cleandata.drop_duplicates(keep='first', inplace=True)
    cleandata.to_csv(os.path.join(datapath,'proteins_no_wiki.tsv'),sep = '\t', header=True)
    print(len(cleandata))
    
#### Generate table of all human genes with a Wikipedia article
def get_genes_with_wiki(datapath):
    url = 'https://query.wikidata.org/sparql'
    query = """
        SELECT ?item ?itemLabel ?geneID ?proteinwdid ?sitelink
        WHERE
        {
          ?item wdt:P31 wd:Q7187 .
          ?item wdt:P703 wd:Q15978631 .
          ?item wdt:P351 ?geneID .
          ?item wdt:P688 ?proteinwdid .
          ?sitelink schema:about ?item .
          FILTER EXISTS {
            ?article schema:about ?item .
            ?article schema:isPartOf <https://en.wikipedia.org/> .
          }

        SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en" }
        }
    """
    r = httprequests.get(url, params = {'format': 'json', 'query': query})
    data = r.json()
    datadf = pd.DataFrame(data['results']['bindings'])
    datadf['uri'] = [x['value'] for x in datadf['item']]
    datadf['label'] = [x['value'] for x in datadf['itemLabel']]
    datadf['geneID'] = [x['value'] for x in datadf['geneID']]
    datadf['QID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['uri']]
    datadf['proteinuri'] = [x['value'] for x in datadf['proteinwdid']]
    datadf['proteinID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['proteinuri']]
    datadf['wikilink'] = [x['value'] for x in datadf['sitelink']]
    cleandata = datadf[['QID','label','geneID','proteinID','wikilink']].copy()
    en_only = cleandata.loc[cleandata['wikilink'].str.contains('en.wiki')].copy()
    en_only.drop_duplicates(keep='first',inplace=True)
    en_only.to_csv(os.path.join(datapath,'genes_en_wiki.tsv'),sep = '\t', header=True)
    print(len(en_only))

#### Generate table of all human genes, where the protein encoded by those genes does have a Wikipedia article
def get_genes_proteins_with_wiki(datapath):
    url = 'https://query.wikidata.org/sparql'
    query = """
        SELECT ?item ?itemLabel ?geneID ?proteinwdid ?sitelink
        WHERE
        {
          ?item wdt:P31 wd:Q7187 .
          ?item wdt:P703 wd:Q15978631 .
          ?item wdt:P351 ?geneID .
          ?item wdt:P688 ?proteinwdid .
          ?sitelink schema:about ?proteinwdid .
          FILTER EXISTS {
            ?article schema:about ?proteinwdid .
            ?article schema:isPartOf <https://en.wikipedia.org/> .
          }

        SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en" }
        }
    """
    r = httprequests.get(url, params = {'format': 'json', 'query': query})
    data = r.json()
    datadf = pd.DataFrame(data['results']['bindings'])
    datadf['uri'] = [x['value'] for x in datadf['item']]
    datadf['label'] = [x['value'] for x in datadf['itemLabel']]
    datadf['geneID'] = [x['value'] for x in datadf['geneID']]
    datadf['QID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['uri']]
    datadf['proteinuri'] = [x['value'] for x in datadf['proteinwdid']]
    datadf['proteinID'] = [x.replace('http://www.wikidata.org/entity/','') for x in datadf['proteinuri']]
    datadf['wikilink'] = [x['value'] for x in datadf['sitelink']]
    cleandata = datadf[['QID','label','geneID','proteinID','wikilink']].copy()
    en_only = cleandata.loc[cleandata['wikilink'].str.contains('en.wiki')].copy()
    en_only.drop_duplicates(keep='first',inplace=True)
    en_only.to_csv(os.path.join(datapath,'proteins_en_wiki.tsv'),sep = '\t', header=True)
    print(len(en_only))
    
def get_wd_info(datapath):
    print("fetching info from Wikidata: ",datetime.now())
    get_genes_no_wiki(datapath)
    get_gene_proteins_no_wiki(datapath)
    get_genes_with_wiki(datapath)
    get_genes_proteins_with_wiki(datapath)
    print("fetching complete: ",datetime.now())
    
#### Merge table of genes and proteins to identify Genes which do NOT have wikipedia articles, 
#### which encode proteins that do NOT have Wikipedia articles
#### ie - Identify genes which show up on both lists (no gene article, no protein article)
def filter_no_wikis(datapath,resultpath):
    genes_no_wiki = pd.read_csv(os.path.join(datapath,'genes_no_wiki.tsv'),delimiter = '\t', header=0, index_col=0)
    proteins_no_wiki = pd.read_csv(os.path.join(datapath,'proteins_no_wiki.tsv'),delimiter = '\t', header=0, index_col=0)
    no_wiki_merge = pd.concat((genes_no_wiki,proteins_no_wiki),ignore_index=True)
    frequency = no_wiki_merge.groupby('geneID').size().reset_index(name='counts')
    no_gene_protein = frequency.loc[frequency['counts']==2]
    no_gene_protein_info = genes_no_wiki.loc[genes_no_wiki['geneID'].isin(no_gene_protein['geneID'].tolist())]
    no_gene_protein_info.to_csv(os.path.join(resultpath,'genes_with_no_gene_protein_wiki.tsv'),sep='t',header=True)

    
#### Merge table of genes and proteins to identify Genes which have a wikipedia article or
#### the protein encoded by those genes has a Wikipedia article
#### ie - Identify all unique wikilinks and their corresponding gene --  this is the merged list
#### Then pull page length for all Wikipedia articles and merged list
#### Filter out Wikipedia articles which are greater than 10,000 characters in length

def filter_wikis(datapath,resultpath):
    genes_en_wiki = pd.read_csv(os.path.join(datapath,'genes_en_wiki.tsv'),delimiter = '\t', header=0, index_col=0)
    proteins_en_wiki = pd.read_csv(os.path.join(datapath,'proteins_en_wiki.tsv'),delimiter = '\t', header=0, index_col=0)
    en_wiki_merge = pd.concat((genes_en_wiki,proteins_en_wiki),ignore_index=True)
    unique_wikis = en_wiki_merge.groupby(['geneID','proteinID','wikilink']).size().reset_index(name='counts')
    unique_wikis.to_csv(os.path.join(datapath,'gene_protein_wikilinks.tsv'),sep='\t',header=True)
    unique_wikis = pd.read_csv(os.path.join(datapath,'gene_protein_wikilinks.tsv'),delimiter='\t',header=0,index_col=0)
    unique_wikis['title'] = [x.replace(" ","_").replace("https://","http://").replace("http://en.wikipedia.org/wiki/","") for x in unique_wikis['wikilink']]
    titlelist = unique_wikis['title'].unique().tolist()
    pageinfo,pagefails = get_wiki_volume_info(mwsite,titlelist)
    wikiinfo = pd.DataFrame(pageinfo)
    wikiinfo.to_csv(os.path.join(datapath,'gene_wiki_vol_info.tsv'),sep='\t',header=True)
    shorter_articles = wikiinfo.loc[wikiinfo['page_length']<10000].copy()
    shorter_articles.sort_values('page_length',ascending=True,inplace=True)
    detailed_shorter_articles = shorter_articles.merge(unique_wikis,on='title',how='inner')
    detailed_shorter_articles.to_csv(os.path.join(resultpath,'priority_by_size.tsv'),sep='\t',header=True)
    print(len(detailed_shorter_articles))
    print(detailed_shorter_articles.head(n=2))
    

%%time
#### Paths
script_path = pathlib.Path(__file__).parent.absolute()
datapath = os.path.join(script_path,'data/')
resultpath = os.path.join(script_path,'results/')

#### Config
useragent = os.environ['User-Agent']

mwsite = mw.Site('en.wikipedia.org', clients_useragent=useragent['User-Agent'])

print("script started: ",datetime.now())
get_wd_info(datapath)
print("wdinfo_script_complete: ",datetime.now())
filter_no_wikis(datapath,resultpath)
print("genes/proteins, no wiki, filtered: ",datetime.now())
filter_wikis(datapath,resultpath)
print("scripts completely run: ",datetime.now())