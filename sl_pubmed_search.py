import streamlit as st
import requests
try:
    from Bio import Entrez
    from Bio import Medline
except ImportError:
    st.error("Biopython is not installed. Please install it using 'pip install biopython'")
    st.stop()

# Set your email for NCBI's E-utilities
Entrez.email = "csv610@gmail.com"

def search_pubmed(query, max_results=10):
    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        # Fetch details for each result
        id_list = record["IdList"]
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        results = []
        for record in records:
            title = record.get("TI", "No title available")
            authors = ", ".join(record.get("AU", ["No authors listed"]))
            journal = record.get("JT", "No journal listed")
            pubdate = record.get("DP", "No date available")
            pmid = record.get("PMID", "No PMID available")
            abstract = record.get("AB", "No abstract available")
            
            results.append({
                "title": title,
                "authors": authors,
                "journal": journal,
                "pubdate": pubdate,
                "pmid": pmid,
                "abstract": abstract
            })
        
        return results
    except Exception as e:
        st.error(f"An error occurred while searching PubMed: {str(e)}")
        return []

# Streamlit app
st.sidebar.title("PubMed")
max_results = st.sidebar.slider("Maximum number of results", 1, 100, 10)

# User input
query = st.text_input("Enter your search query:")

if query:  # Trigger search when user presses "Return"
    with st.spinner("Searching PubMed..."):
        results = search_pubmed(query, max_results)
    
    if results:
        for result in results:
            st.write(f"**Title:** {result['title']}")
            st.write(f"**Authors:** {result['authors']}")
            st.write(f"**Journal:** {result['journal']}")
            st.write(f"**Publication Date:** {result['pubdate']}")
            st.write(f"**PMID:** {result['pmid']}")
            with st.expander("Show Abstract"):
                st.write(result['abstract'])
            st.write("---")
    else:
        st.warning("No results found or an error occurred. Please try again.")
