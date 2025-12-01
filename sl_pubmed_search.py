import streamlit as st
from search_engines import PubMedSearch


@st.cache_resource
def get_search_engine():
    return PubMedSearch()


st.sidebar.title("PubMed")
max_results = st.sidebar.slider("Maximum number of results", 1, 100, 10)

query = st.text_input("Enter your search query:")

if query:
    with st.spinner("Searching PubMed..."):
        try:
            search_engine = get_search_engine()
            results = search_engine.search(query, max_results)

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
        except Exception as e:
            st.error(f"Error: {e}")
