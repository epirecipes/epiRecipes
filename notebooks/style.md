## Style guide for notebooks

1. Make a separate folder for each model.
2. Background for each model should be placed into a notebook `intro.ipynb`, which will then be converted into a Markdown document, `intro.md`.
3. Implementations for each model should go into a separate file:
  - `r.ipynb`: R
  - `julia.ipynb`: Julia
  - `python.ipynb`: Python
4. If there are separate implementations for each language, use a suffix separated by a hyphen e.g. `r-simecol.ipynb`.
5. Each implementation should follow the following structure:
   - Title of implementation
   - Name of author followed by GitHub username
   - Date in IS)8601
   - Simulation
   - Visualisation
6. Ensure there are no empty cells at the end of the notebook.
7. Other style points:
   - Use `snake_case` for variables.
   - Use `##` for headings and `###` for subheadings.
