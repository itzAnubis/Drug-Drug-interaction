import { ApiResponse, Interaction } from '../types/interactions';

// Update this URL to your actual API endpoint
const API_URL = "API_URL"; // put here the API from the model_api_DDI.ipynb

export const predictInteractions = async (
  drug1Name: string,
  drug2Name: string
): Promise<{
  interactions: Interaction[];
  drug1_smiles: string;
  drug2_smiles: string;
  drug1_img_base64: string;
  drug2_img_base64: string;
}> => {
  try {
    const response = await fetch(API_URL, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json'
      },
      body: JSON.stringify({
        drug1: drug1Name,
        drug2: drug2Name,
      }),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(
        `API request failed: ${response.status} ${response.statusText}${
          errorText ? ` - ${errorText}` : ''
        }`
      );
    }

    const data: ApiResponse = await response.json();

    // Convert the API response to our internal format
    const interactions: Interaction[] = [{
      type: 'Interaction',
      severity: determineSeverity(data.interaction),
      description: data.interaction,
    }];

    return {
      interactions,
      drug1_smiles: data.drug1_smiles,
      drug2_smiles: data.drug2_smiles,
      drug1_img_base64: data.drug1_img_base64,
      drug2_img_base64: data.drug2_img_base64,
    };
  } catch (error) {
    console.error('Error predicting interactions:', error);
    if (error instanceof Error) {
      // Provide more specific error message based on the type of error
      if (error.message.includes('Failed to fetch')) {
        throw new Error('Unable to connect to the API server. Please check if the server is running and try again.');
      }
      throw new Error(`API Error: ${error.message}`);
    }
    throw error;
  }
};

// Simple function to determine severity based on keywords in the interaction description
function determineSeverity(description: string): 'High' | 'Moderate' | 'Low' {
  const lowercaseDesc = description.toLowerCase();
  if (
    lowercaseDesc.includes('severe') ||
    lowercaseDesc.includes('dangerous') ||
    lowercaseDesc.includes('fatal') ||
    lowercaseDesc.includes('critical')
  ) {
    return 'High';
  } else if (
    lowercaseDesc.includes('moderate') ||
    lowercaseDesc.includes('significant') ||
    lowercaseDesc.includes('notable')
  ) {
    return 'Moderate';
  }
  return 'Low';
}