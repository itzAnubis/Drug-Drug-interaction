import nest_asyncio
import uvicorn
import threading
from pyngrok import ngrok
import requests

nest_asyncio.apply()

try:
    public_url = ngrok.connect(8000)
    print(f"🚀 Your public API is available at: {public_url}/docs")

    sheet_url = "https://script.google.com/macros/s/AKfycbw14Qm1u1BlV9fZp83h0_7GGKLTDk9arMx6MHvHw4jZrt7A2zzBBS0yLdawINNqPM2R/exec"
    requests.post(sheet_url, json={"url": str(public_url)})

except Exception as e:
    print("❌ Failed to start Ngrok tunnel:", str(e))

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from transformers import T5ForConditionalGeneration, T5Config, T5Tokenizer
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp
import torch
import base64
from io import BytesIO

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
config = T5Config.from_pretrained("t5-small")
model = T5ForConditionalGeneration(config)

try:
    model.load_state_dict(torch.load("/content/drive/MyDrive/Copy of t5_interaction_gen.pt", map_location=device))
    print("✅ Model loaded successfully")
except Exception as e:
    print("❌ Error loading model:", str(e))

tokenizer = T5Tokenizer.from_pretrained("t5-small")
model.to(device)
model.eval()

class DrugPair(BaseModel):
    drug1: str
    drug2: str

def name_to_smiles(query: str):
    mol = Chem.MolFromSmiles(query)
    if mol is not None:
        return query
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

# Run FastAPI app using threading (for Colab)
def run_app():
    uvicorn.run(app, host="0.0.0.0", port=8000)

thread = threading.Thread(target=run_app)
thread.start()
