import streamlit as st
from search_engines import DuckDuckGoSearch


@st.cache_resource
def get_search_engine():
    return DuckDuckGoSearch()


st.set_page_config(layout="wide", page_title="DuckDuckGo Search")

st.sidebar.title("DuckDuckGo")
num_results = st.sidebar.slider("Number of results:", min_value=1, max_value=100, value=10)
safesearch = st.sidebar.selectbox("SafeSearch level:", ["Off", "Moderate", "Strict"], index=1)

search_query = st.text_input("Enter your search query:")

if search_query:
    st.write(f"Searching for: {search_query} with SafeSearch level: {safesearch}")

    try:
        search_engine = get_search_engine()
        results = search_engine.search(search_query, num_results, safesearch=safesearch)

        if results:
            for i, result in enumerate(results, 1):
                st.subheader(f"Result {i}: {result['title']}")
                st.write(f"**Link:** [{result['href']}]({result['href']})")
                st.write(f"**Summary:** {result['body']}")
                st.divider()
        else:
            st.warning("No results found.")
    except Exception as e:
        st.error(f"Error: {e}")
