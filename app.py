from fastapi import FastAPI
from pydantic import BaseModel
from transformers import T5ForConditionalGeneration, T5Config, T5Tokenizer
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp
import torch
import base64
from io import BytesIO

app = FastAPI()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")



config = T5Config.from_pretrained("t5-small")  # Or custom config path if you have one
model = T5ForConditionalGeneration(config)
model.load_state_dict(torch.load("models/t5_interaction_gen.pt", map_location=device))
tokenizer = T5Tokenizer.from_pretrained("t5-small")
model.to(device)
model.eval()


class DrugPair(BaseModel):
    drug1: str
    drug2: str

def name_to_smiles(query: str):
    """Convert drug name to SMILES using PubChem if needed"""
    mol = Chem.MolFromSmiles(query)
    if mol is not None:
        return query  # It's already a valid SMILES
    try:
        compound = pcp.get_compounds(query, 'name')[0]
        return compound.isomeric_smiles
    except Exception:
        return None

def mol_to_base64(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=(224, 224))
    buf = BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("utf-8")

@app.post("/predict")
def predict_interaction(data: DrugPair):
    drug1 = data.drug1
    drug2 = data.drug2
    
    smi1 = name_to_smiles(drug1)
    smi2 = name_to_smiles(drug2)

    if smi1 is None or smi2 is None:
        return {"error": "One or both drug names could not be converted to SMILES."}

    input_text = f"Drug1: {drug1} Drug2: {drug2}"
    input_ids = tokenizer.encode(input_text, return_tensors="pt").to(device)
    output = model.generate(input_ids, max_length=128)
    result = tokenizer.decode(output[0], skip_special_tokens=True)

    img1 = mol_to_base64(smi1)
    img2 = mol_to_base64(smi2)

    return {
        "interaction": result,
        "drug1_smiles": smi1,
        "drug2_smiles": smi2,
        "drug1_img_base64": img1,
        "drug2_img_base64": img2
    }
