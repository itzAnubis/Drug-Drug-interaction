export interface Drug {
  smiles: string;
  name: string;
  img_base64?: string;
}

export interface Interaction {
  type: string;
  severity: 'High' | 'Moderate' | 'Low';
  description: string;
  recommendation?: string;
}

export interface InteractionResult {
  drug1: Drug;
  drug2: Drug;
  interactions: Interaction[];
  timestamp: number;
}

export interface ApiResponse {
  interaction: string;
  drug1_smiles: string;
  drug2_smiles: string;
  drug1_img_base64: string;
  drug2_img_base64: string;
}