'''
Runs server for the protein landscape computation.
'''

from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import parse_qs
from epicore_utils.modules.visualize_protein import plot_protein_landscape
import pandas as pd
import ast
import time
import numpy as np


class MyHTTPRequestHandler(BaseHTTPRequestHandler):
    ''' HTTPRequestHandler
    '''

    def do_GET(self):
        ''' GET method of the HTTPRequestHandler.
        '''        
        report = getattr(self.server, "report", "No report available")
        self.send_response(200,message='well')
        self.end_headers()
             
        if self.path == '/prot_lan.svg':
            with open('prot_lan.svg', 'rb') as f:
                self.wfile.write(f.read())
        elif 'svg' in self.path:
            with open(self.path[1:], 'rb') as f:
                self.wfile.write(f.read())
        else:
            self.wfile.write(report.encode("utf-8"))


    def do_POST(self):
        ''' POST method of the HTTPRequestHandler.
        '''
        # get accession, epicore result and proteome dictionary
        content_length = int(self.headers.get('Content-Length', 0))
        post_data = self.rfile.read(content_length) 
        parsed_data = parse_qs(post_data.decode('utf-8'))
        accession = parsed_data.get('accession', [''])[0]
        epicore_csv = parsed_data.get('epicore_csv', [''])[0]
        proteome_dict = parsed_data.get('proteome_dict', [''])[0]

        proteome_dict = ast.literal_eval(proteome_dict)

        # read in epicore result
        protein_df = pd.read_csv(epicore_csv)
    
        protein_df['grouped_peptides_start'] = protein_df['grouped_peptides_start'].apply(ast.literal_eval)
        protein_df['core_epitopes_start'] = protein_df['core_epitopes_start'].apply(lambda cell: eval(cell, {"np": np}))
        protein_df['core_epitopes_end'] = protein_df['core_epitopes_end'].apply(lambda cell: eval(cell, {"np": np}))
        protein_df['landscape'] = protein_df['landscape'].apply(ast.literal_eval)

        img_url = f'prot_lan{str(int(time.time()))}.svg'
        # compute protein landscape
        fig = plot_protein_landscape(protein_df, accession, proteome_dict)
        fig.savefig(img_url,bbox_inches='tight')
        
        self.send_response(200,message=accession)
        self.end_headers()
        self.wfile.write(img_url.encode("utf-8"))
        

    def do_OPTIONS(self):
        ''' OPTIONS method of the HTTPRequestHandler.
        '''
        self.send_response(200, "ok")
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET,POST, OPTIONS')
        self.send_header("Access-Control-Allow-Headers", "X-Requested-With")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()


def run(port,report, server_class=HTTPServer,handler_class=MyHTTPRequestHandler):
    ''' Runs the server.

    Args:
        port: An integer defining the port.
        server_class: Specifies which server class is used.
        handler_class: Specifies which request handler class is used.
    
    '''
    server_address = ('', port)
    httpd = server_class(server_address, handler_class)
    httpd.report = report 
    httpd.serve_forever()