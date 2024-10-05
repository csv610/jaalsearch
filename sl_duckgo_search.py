import streamlit as st
from duckduckgo_search import DDGS

# Cache the DDGS instance to avoid creating multiple instances
@st.cache_resource
def get_ddgs_instance():
    return DDGS()

# Perform a search using DuckDuckGo (removed caching due to unhashable parameter)
def search_duckduckgo(ddgs, query, num_results=5, safesearch="Moderate"):
    try:
        return list(ddgs.text(query, max_results=num_results, safesearch=safesearch))
    except Exception as e:
        st.error(f"An error occurred while searching: {e}")
        return []

# Callback function for handling input changes
def on_search_query_change():
    st.session_state['search'] = True

# Set the Streamlit page configuration
st.set_page_config(layout="wide", page_title="DuckDuckGo Search")

# Get the DDGS instance
ddgs_instance = get_ddgs_instance()

# Sidebar for search settings
st.sidebar.title("DuckDuckGo")
num_results = st.sidebar.slider("Number of results:", min_value=1, max_value=100, value=10)
safesearch = st.sidebar.selectbox("SafeSearch level:", ["Off", "Moderate", "Strict"], index=1)

# Main page title

# Input for search query
search_query = st.text_input("Enter your search query:", key="search_query", on_change=on_search_query_change)

# Perform search and display results
if search_query:
    st.write(f"Searching for: {search_query} with SafeSearch level: {safesearch}")
    
    results = search_duckduckgo(ddgs_instance, search_query, num_results, safesearch=safesearch)
    
    if results:
        for i, result in enumerate(results, 1):
            st.subheader(f"Result {i}: {result['title']}")
            st.write(f"**Link:** [{result['href']}]({result['href']})")
            st.write(f"**Summary:** {result['body']}")
            st.divider()
    else:
        st.warning("No results found.")
