import streamlit as st
from search_engines import WikipediaSearch


@st.cache_resource
def get_search_engine():
    return WikipediaSearch()


st.sidebar.title("Wikipedia Search")

num_results = st.sidebar.slider("Number of results", 1, 50, 10)

query = st.text_input("Enter your search query and press Enter:")

if query:
    with st.spinner("Searching Wikipedia..."):
        try:
            search_engine = get_search_engine()
            results = search_engine.search(query, num_results)

            if results:
                for result in results:
                    st.markdown(f"### {result['title']}")
                    if result['url']:
                        st.markdown(f"[Link to Wikipedia page]({result['url']})")
                    st.write(result['summary'])
                    st.divider()
            else:
                st.warning(f"No results found for '{query}'.")
        except Exception as e:
            st.error(f"Error: {e}")
