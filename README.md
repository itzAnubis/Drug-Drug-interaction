
# Drug-Drug Interaction (DDI) Prediction System

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-orange)](https://pytorch.org/)
[![Transformers](https://img.shields.io/badge/Transformers-T5--small-yellow)](https://huggingface.co/docs/transformers)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.95%2B-green)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/React-18%2B-blue)](https://react.dev/)

```markdown
A machine learning system that predicts potential interactions between two drugs using a fine-tuned T5-small language model, complete with molecular visualization.

## Table of Contents
- [Features](#features)
- [Demo](#demo)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Model Details](#model-details)
- [API Endpoints](#api-endpoints)
- [Evaluation](#evaluation)
- [Future Work](#future-work)
- [Limitations](#limitations)
- [Contributing](#contributing)

## Features
- **Drug Interaction Prediction**: Predicts potential interactions between any two drugs
- **Molecular Visualization**:
  - Displays 2D molecular structures
  - Shows SMILES notations
- **REST API**: FastAPI backend with Swagger documentation
- **Web Interface**: React-based frontend for easy interaction
- **Chemical Validation**: Uses PubChemPy and RDKit for chemical structure handling

## Demo
**Sample Input:**
```plaintext
Drug 1: Aspirin
Drug 2: Paracetamol
```

**Sample Output:**
```plaintext
Predicted Interaction: The risk or severity of adverse effects can be increased when Paracetamol is combined with Aspirin.

SMILES Notations:
- Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
- Paracetamol: CC(=O)NC1=CC=C(C=C1)O
```

## Installation

### Backend Setup
```bash
# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Frontend Setup
```bash
cd frontend
npm install
npm run dev
```

### Requirements
Backend requirements (`requirements.txt`):
```
fastapi==0.95.2
uvicorn==0.22.0
torch==2.0.1
transformers==4.30.2
pubchempy==1.0.4
rdkit==2023.9.5
numpy==1.26.4
sentencepiece==0.1.99
python-multipart==0.0.6
```

## Usage

### Running the Backend
**Option 1:** Using Python script
```bash
uvicorn main:app --reload
```
API will be available at `http://localhost:8000`

**Option 2:** Using Jupyter Notebook
1. Open `model_api_DDI.ipynb`
2. Run all cells to start the API
3. Update the frontend's API URL in `frontend/src/services/interactionService.ts`

### Running the Frontend
```bash
cd frontend
npm start
```
Frontend will launch at `http://localhost:5173`

### Using the API Directly
```bash
curl -X POST "http://localhost:8000/predict-interaction" \
-H "Content-Type: application/json" \
-d '{"drug1": "Aspirin", "drug2": "Paracetamol"}'
```

## Project Structure
```
drug-drug-interaction/
├── frontend/                  # React frontend application
│   ├── public/                # Static assets
│   └── src/                   # React source files
├── model/                     # Model files
│   ├── t5-small/              # Pretrained model weights
├── drug-drug-interactions.csv # Training dataset (1,634 drug pairs)
├── model_api_DDI.ipynb        # Model development notebook
├── main.py                    # FastAPI application
├── requirements.txt           # Python dependencies
└── README.md                  # This file
```

## Model Details
- **Base Model**: T5-small (Text-to-Text Transfer Transformer)
- **Input Format**: "drug1: {drug1_name}; drug2: {drug2_name}"
- **Output Format**: Natural language interaction description
- **Training Data**: 1,634 drug pairs from drug-drug-interactions.csv
- **Training Hardware**: NVIDIA T4 GPU (Google Colab)
- **Tokenization**: SentencePiece tokenizer

## API Endpoints
- `POST /predict-interaction` - Predict interaction between two drugs
  - Request body: `{"drug1": string, "drug2": string}`
  - Response: 
    ```json
    {
      "interaction": string,
      "drug1_smiles": string,
      "drug2_smiles": string,
      "drug1_img": "base64_encoded_image",
      "drug2_img": "base64_encoded_image"
    }
    ```
- `GET /docs` - Interactive API documentation (Swagger UI)

## Evaluation
The model was evaluated using:
- **Loss Function**: Cross-entropy loss
- **Training Loss**: 0.015900
- **Validation Strategy**: 80-20 train-test split

## Future Work
- [ ] Expand drug coverage to 10,000+ medications
- [ ] Add interaction severity scoring (mild/moderate/severe)
- [ ] Implement batch prediction for multiple drug pairs
- [ ] Add user accounts and interaction history
- [ ] Develop mobile application version
- [ ] Incorporate drug dosage considerations
- [ ] Add literature references for interactions

## Limitations
1. Currently covers ~1,800 common medications
2. Accuracy may vary for rare drug combinations
3. Doesn't account for patient-specific factors (age, weight, etc.)
4. Predictions should be verified by healthcare professionals
5. Limited to pairwise interactions (not accounting for multiple drug combinations)

## Contributing
Contributions are welcome! Please follow these steps:
1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

**Disclaimer**: This system provides predictions only and should not be used as a substitute for professional medical advice. Always consult a healthcare provider about potential drug interactions.

