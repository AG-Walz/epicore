'''
Generates a report for the results of epicore.
'''

import webbrowser 
import ast
import pandas as pd
import threading
import epicore_utils.modules.locserver


def gen_report(length_distribution: str, intensity_hist: str, epitope_df: pd.DataFrame, peps: int, epitopes: int, n_removed_peps: int, ctx, evidence_file: str, epicore_csv: str):
    """ Generates a report including the most important information of the results. 

    Args:
        length_distribution: The path to the file containing the length 
            distribution of peptides and epitopes. 
        intensity_hist: The path to the file containing a histogram of the 
            number of peptides contributing to each core. 
        epitope_df: A dataframe containing one epitope and information of that  
            epitope per row. 
        peps: Number of peptides in the evidence file. 
        epitopes: Number of epitopes generated.
        n_removed_peps: Number of peptides removed from the evidence file due 
            to the absence of their accessions in the proteome.
        ctx: Object containing the input parameters.
        evidence_file: Path to the evidence file.
        epicore_csv: Path to the epicore result csv.

    Returns:
        Open a html report summarizing some information of the epicore 
        result. For example the report includes a figure of the peptide/epitope 
        length distribution, a histogram of the number of peptides mapped to 
        each epitope and the ten epitopes with the highest number of mapped 
        peptides. 
    """

    # compute mean of epitope core sequence
    len_core_epitopes = epitope_df['consensus_epitopes'].apply(lambda sequence: len(sequence))
    mean_core_length = sum(len_core_epitopes) / len(len_core_epitopes)

    # compute mean of epitope whole sequence
    len_whole_epitopes = epitope_df['whole_epitopes'].apply(lambda sequence: len(sequence))
    mean_whole_length = sum(len_whole_epitopes) / len(len_whole_epitopes)

    epitope_df_sort = epitope_df.sort_values(by='landscape', key=lambda landscapes: landscapes.apply(lambda landscape: max(ast.literal_eval(landscape))), ascending=False).head(10)
    epitope_df_sort_html = epitope_df_sort.to_html()

    # get accession of top ten epitopes to add them as options to the combobox
    epitope_df_sort['accession'] = epitope_df_sort['accession'].apply(lambda row: row.split(','))
    epitope_df_sort_accessions = epitope_df_sort['accession'].explode().tolist() 

    html_content = f''' <!DOCTYPE html>
                        <html>
                        <head>
                        <script src="https://cdn.jsdelivr.net/npm/brython@3/brython.min.js">
                        </script>
                        <script src="https://cdn.jsdelivr.net/npm/brython@3/brython_stdlib.js">
                        </script>
                            <style> 
                                .row {{
                                    display: flex;
                                    justify-content: space-between;
                                }}
                                .column {{
                                    flex: 1;
                                    margin: 0 10px;
                                }}
                                h1{{text-align: center;}}
                            </style>
                            <title>Epicore report</title>
                        </head>
                        <body>
                            <h1>Epicore report</h1>
                            <div class="row">
                                <div class="column"><p> The histogram shows the number of peptides/epitopes of a certain length. {peps} peptides were reduced to {epitopes} epitopes.</p><img src="{length_distribution}"></div>
                                <div class="column"><p> The histogram visualizes how many peptides contribute to the different epitopes.<br></p><img src="{intensity_hist}"></div>
                            </div>

                            <script type="text/python">
                                from browser import document, bind, ajax, html

                                select_element = document["exampleList"]
                                for item in {epitope_df_sort_accessions}:
                                    option = html.OPTION(value=item)
                                    select_element <= option

                                def on_complete(req):
                                    if req.status == 200 or req.status == 0:
                                        image_url = req.text  
                                        img_element = document["image_container"]
                                        img_element.clear()
                                        img_element <= html.IMG(src=image_url, alt="The accession you selected is not in your input data. Try a different accession.", style="width: 1400px;")

                                @bind('#subm','click')
                                def welcome(event):
                                    event.preventDefault()
                                    req = ajax.Ajax()
                                    req.bind('complete', on_complete)
                                    accession = document['fname'].value
                                    req.open('POST','http://localhost:8000/', True)
                                    req.send({{'accession': accession, 'epicore_csv':'{epicore_csv}','proteome_dict':{ctx.obj.proteome_dict}}})
 
                            </script>

                            <form id="prot_landscape">
                                <label for="fname">Protein accession: </label>
                                <input type="text" id="fname" name="fname" list="exampleList"><br>
                                <datalist id="exampleList">
                                </datalist>
                                <button id="subm">Plot protein landscape</button>
                            </form> 
                            <div id="image_container"></div>
                            <p id="output"></p>
                            <p> The average length of an epitope core sequence is {mean_core_length}. <br> The average length of an epitope sequence is {mean_whole_length}.</p><br>
                            <p>{n_removed_peps} peptides were removed, since they do not appear in the provided evidence file.</p>
                            <p> The used parameter values:<br> - min_epi_length:{ctx.obj.min_epi_length}<br> 
                            - min_overlap:{ctx.obj.min_overlap}<br> - max_step_size:{ctx.obj.max_step_size}<br> - seq_column:{ctx.obj.seq_column}<br>
                            - protacc_column:{ctx.obj.protacc_column} <br> - intensity_column:{ctx.obj.intensity_column} <br> - start_column: {ctx.obj.start_column} <br>
                            - end_column: {ctx.obj.end_column} <br> - out_dir: {ctx.obj.out_dir} <br> - mod_pattern: {ctx.obj.mod_pattern} <br> - delimiter: {ctx.obj.delimiter} <br>
                            The input evidence file: {evidence_file}</p>
                            <p> The ten epitopes with the highest number of mapped peptides:<br>{epitope_df_sort_html}</p>
                            

                        </body>
                        </html>'''
    
    threading.Thread(target=epicore_utils.modules.locserver.run, args=(8000,html_content)).start()
    webbrowser.open_new_tab('http://localhost:8000/report.html')