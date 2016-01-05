#!/usr/bin/python
import json
from flask import Flask,request,jsonify,make_response
from flask_restful import Resource, Api
from sqlalchemy import create_engine

e = create_engine('mysql+pymysql://root:paok@localhost/xwiki')

app = Flask(__name__)
api = Api(app)

class HPO_Data(Resource):
    def post(self,hpo):
         
        conn = e.connect()

        hp_list=[]
        ids_list =[]
        fin_ids = {}

        #Get JSON from POST method
        hp_json = request.json   

        #From JSON get all the "features"(HPO) objects
        features = hp_json['features']

        #Append all HPO terms to a list
        for hp_dict in features:
        	for key in hp_dict:
        		hp_list.append(hp_dict[key])
        
        #For each HPO term get the IDs from the database and store them in list
        for i in hp_list:
            query1 = conn.execute("SELECT GROUP_CONCAT(XWL_ID separator ', ') FROM XWIKILISTITEMS WHERE XWL_NAME IN ('phenotype','extended_phenotype') AND XWL_VALUE='%s'"%i)
            ids_list.append(list(query1.cursor.fetchall()))
       
        #Create a dictionary that as keys, has the HPO terms and as values lists of all the corresponding IDs
        for j in range(len(ids_list)):
        	fin_ids[hp_list[j]] = ids_list[j]


        #Helping data structures for next part
        list1 = []          
        r={}

       
        #For every HPO term in our dictionary we query the database with the ID values in order to get the corresponding "Sample" name from another table (ID is the key,same between the tables)
        for key in fin_ids:
            cleaned=[i[0] for i in fin_ids[key]]    
            for j in cleaned:
              #We query only if we have an ID, meaning that the actual HPO term posted in JSON,exists in the database.
              if j is not None:
                query3 = conn.execute("SELECT XWS_VALUE as external_id FROM XWIKISTRINGS WHERE XWS_NAME = 'external_id' AND XWS_ID IN (" + str(j) + ")")
                external = query3.cursor.fetchone()
                while external is not None:
                  external=str(external).translate(None, "(),'")
                  external=external.replace(" ", ',')
                  list1.append({external : key})
                  external=query3.cursor.fetchone()
      
        for d in list1:
          ((x, y),) = d.items() 
          r.setdefault(x, []).append(y)

        vam = [ {k: v} for (k, v) in r.items() ]
        
        #Dict with patient data
        patients = {'patient': vam }
        
        #Final dict with all data
        result = {'results' : patients}

        return result
        


class Patient_Data(Resource):
    def get(self, ext_id):

        conn = e.connect()

        #We either parse the query parameter 'eid' and feed it to the query or directly we put the ext_id
    	ext_id2=request.args.get("id")


        query1 = conn.execute("SELECT XWS_ID FROM XWIKISTRINGS WHERE XWS_VALUE='%s'"%ext_id2)     
        xws_id = str(query1.fetchall())
        xws_id = xws_id.translate(None, '()[],')

        if (xws_id != ''):
   
              name = conn.execute("SELECT (SELECT XWS_VALUE FROM XWIKISTRINGS WHERE XWS_NAME = 'first_name' AND XWS_ID=" + xws_id + ") as first_name, + \
              (SELECT XWS_VALUE FROM XWIKISTRINGS WHERE XWS_NAME = 'last_name' AND XWS_ID=" + xws_id + ") as last_name")       
              

              ethnicity = conn.execute("SELECT (SELECT GROUP_CONCAT(XWL_VALUE separator ', ') from XWIKILISTITEMS WHERE XWL_NAME='maternal_ethnicity' AND XWL_ID =" + xws_id + ") as maternal_ethnicity, + \
              (SELECT GROUP_CONCAT(XWL_VALUE separator ', ') from XWIKILISTITEMS WHERE XWL_NAME='paternal_ethnicity' AND XWL_ID =" + xws_id + ") as paternal_ethnicity")
  

              external_id= conn.execute("SELECT XWS_VALUE as ex_id "+ \
                                   "FROM XWIKISTRINGS " + \
                                   "WHERE XWS_ID =" + xws_id + " AND XWS_NAME = 'external_id'")

              gender = conn.execute("SELECT XWS_VALUE as gender "+ \
                                   "FROM XWIKISTRINGS " + \
                                   "WHERE XWS_ID =" + xws_id + " AND XWS_NAME = 'gender'")

              phenotype = conn.execute("SELECT DISTINCT  XWL_VALUE as id "+ \
                                   "FROM XWIKILISTITEMS " + \
                                   "WHERE XWL_ID =" + xws_id + " AND XWL_NAME = 'phenotype'")


              extended_phenotype = conn.execute("SELECT DISTINCT  XWL_VALUE as id "+ \
                                   "FROM XWIKILISTITEMS " + \
                                   "WHERE XWL_ID =" + xws_id + " AND XWL_NAME = 'extended_phenotype'")


              prenatal_phenotype = conn.execute("SELECT DISTINCT  XWL_VALUE as id "+ \
                                   "FROM XWIKILISTITEMS " + \
                                   "WHERE XWL_ID =" + xws_id + " AND XWL_NAME = 'prenatal_phenotype'")
        
              
              disorders = conn.execute("SELECT XWL_VALUE as OMIM_id "+ \
                                   "FROM XWIKILISTITEMS " + \
                                   "WHERE XWL_ID =" + xws_id + " AND XWL_NAME = 'omim_id'")


              gene = conn.execute("SELECT XWS_VALUE as gene "+ \
                                   "FROM XWIKISTRINGS " + \
                                   "WHERE XWS_ID =" + xws_id + " AND XWS_NAME = 'gene'")
              
              birth = conn.execute("SELECT XWS_VALUE as b_date "+ \
                                   "FROM XWIKISTRINGS " + \
                                   "WHERE XWS_ID =" + xws_id + " AND XWS_NAME = 'date_of_birth_entered'")

              death = conn.execute("SELECT XWS_VALUE as d_date "+ \
                                   "FROM XWIKISTRINGS " + \
                                   "WHERE XWS_ID =" + xws_id + " AND XWS_NAME = 'date_of_death_entered'")

              age_onset = conn.execute("SELECT XWS_VALUE as id "+ \
                                   "FROM XWIKISTRINGS " + \
                                   "WHERE XWS_ID =" + xws_id + " AND XWS_NAME = 'global_age_of_onset'")

              inheritance = conn.execute("SELECT XWL_VALUE as id "+ \
                                   "FROM XWIKILISTITEMS " + \
                                   "WHERE XWL_ID =" + xws_id + " AND XWL_NAME = 'global_mode_of_inheritance'")

            
              result = {'features': [dict(zip(tuple (phenotype.keys()) ,i)) for i in phenotype.cursor],
                        'extended_features': [dict(zip(tuple (extended_phenotype.keys()) ,i)) for i in extended_phenotype.cursor],
                        'patient_name': [dict(zip(tuple (name.keys()) ,i)) for i in name.cursor],
                        'genes': [dict(zip(tuple (gene.keys()) ,i)) for i in gene.cursor],
                        'global_age_of_onset': [dict(zip(tuple (age_onset.keys()) ,i)) for i in age_onset.cursor],
                        'ethnicity': [dict(zip(tuple (ethnicity.keys()) ,i)) for i in ethnicity.cursor],
                        'disorders': [dict(zip(tuple (disorders.keys()) ,i)) for i in disorders.cursor],
                        'prenatal_perinatal_phenotype': [dict(zip(tuple (prenatal_phenotype.keys()) ,i)) for i in prenatal_phenotype.cursor],
                        'global_mode_of_inheritance': [dict(zip(tuple (inheritance.keys()) ,i)) for i in inheritance.cursor],
                        'sex': str(gender.cursor.fetchone()).translate(None, "'()[],'"),
                        'external_id': str(external_id.cursor.fetchone()).translate(None, "'()[],'"),
                        'date_of_birth': str(birth.cursor.fetchone()).translate(None, "'()[],'"),
                        'date_of_death': str(death.cursor.fetchone()).translate(None, "'()[],'")} 

              
              #Replace all "None" values of dictionary with empty string
              for key in result :
                  if result[key] == "None":
                     result[key] = ""


              #Handling of multiple SELECT commands to manage null values
              ethn_none = 0
              name_none = 0

              result["ethnicity"] = result.get("ethnicity")[0]
              result["patient_name"] = result.get("patient_name")[0]


              for key,val in result["ethnicity"].items():
                 if val is None:
                    result["ethnicity"][key] = ""
                    ethn_none += 1


              for key,val in result["patient_name"].items():
                 if val is None:
                    result["patient_name"][key] = "" 
                    name_none += 1


              if ethn_none > 1:
                 result["ethnicity"] = {}


              if name_none > 1:
                 result["patient_name"] = {}

          
              #Remove keys with empty values from dictionary 
              #result=dict((k, v) for k, v in result.iteritems() if v)
       
                                
        else:
           result=make_response(jsonify({'error': 'Identifier not found'}), 404)
        
        
        return result

api.add_resource(Patient_Data, '/patients/<string:ext_id>')
api.add_resource(HPO_Data, '/HPO/<string:hpo>')

if __name__ == '__main__':
     app.run(debug=True)
