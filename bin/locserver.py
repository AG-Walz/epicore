'''
Runs server for the protein landscape computation.
'''

from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import parse_qs
from bin.visualize_protein import vis_prot
import pandas as pd
import ast



class MyHTTPRequestHandler(BaseHTTPRequestHandler):
    ''' HTTPRequestHandler
    '''

    def do_GET(self):
        ''' GET method of the HTTPRequestHandler.
        '''        
        report = getattr(self.server, "report", "No report available")
        self.send_response(200,message='well')
        self.end_headers()
             
        if self.path == '/prot_lan.png':
            with open('prot_lan.png', 'rb') as f:
                self.wfile.write(f.read())
        elif 'png' in self.path:
            with open(self.path[1:], 'rb') as f:
                self.wfile.write(f.read())
        else:
            self.wfile.write(report.encode("utf-8"))


    def do_POST(self):
        ''' POST method of the HTTPRequestHandler.
        '''
        # get accession, plateau result and proteome dictionary
        content_length = int(self.headers.get('Content-Length', 0))
        post_data = self.rfile.read(content_length) 
        parsed_data = parse_qs(post_data.decode('utf-8'))
        accession = parsed_data.get('accession', [''])[0]
        plateau_csv = parsed_data.get('plateau_csv', [''])[0]
        proteome_dict = parsed_data.get('proteome_dict', [''])[0]

        proteome_dict = ast.literal_eval(proteome_dict)

        # read in plateau_result
        protein_df = pd.read_csv(plateau_csv)
        protein_df['grouped_peptides_start'] = protein_df['grouped_peptides_start'].apply(ast.literal_eval)
        protein_df['core_epitopes_start'] = protein_df['core_epitopes_start'].apply(ast.literal_eval)
        protein_df['core_epitopes_end'] = protein_df['core_epitopes_end'].apply(ast.literal_eval)
        protein_df['landscape'] = protein_df['landscape'].apply(ast.literal_eval)

        # compute protein landscape
        fig = vis_prot(protein_df, accession, proteome_dict)
        fig.savefig('prot_lan.png')
        
        self.send_response(200,message=accession)
        self.end_headers()
        self.wfile.write('/prot_lan.png'.encode("utf-8"))
        

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

