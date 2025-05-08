import { Interaction } from '../types/interactions';

// Mock service to simulate API calls to a drug interaction prediction model
export const mockPredictInteractions = async (
  drug1Smiles: string, 
  drug2Smiles: string
): Promise<Interaction[]> => {
  // Simulate network delay
  await new Promise(resolve => setTimeout(resolve, 1500));
  
  // Map some known drug SMILES to common names for our mock
  const knownDrugs: Record<string, string> = {
    'CC(=O)OC1=CC=CC=C1C(=O)O': 'Aspirin',
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O': 'Ibuprofen',
    'CC(=O)NC1=CC=C(C=C1)O': 'Paracetamol',
    'CC(=O)CC(C1=CC=CC=C1)C1=C(O)C2=CC=CC=C2OC1=O': 'Warfarin',
    'CN1C=NC2=C1C(=O)N(C)C(=O)N2C': 'Caffeine',
    'CN1C2=C(C(=O)N(C)C1=O)N(C)C=N2': 'Theophylline'
  };
  
  // Predefined mock interactions for specific drug pairs
  const knownInteractions: Record<string, Interaction[]> = {
    // Aspirin + Warfarin
    'Aspirin-Warfarin': [
      {
        type: 'Pharmacodynamic',
        severity: 'High',
        description: 'Increased risk of bleeding due to combined anticoagulant effects.',
        recommendation: 'Avoid concurrent use. If combination is necessary, monitor bleeding times closely and consider dose adjustment.'
      },
      {
        type: 'Pharmacokinetic',
        severity: 'Moderate',
        description: 'Aspirin may displace warfarin from protein binding sites, increasing free warfarin levels.',
        recommendation: 'Monitor INR more frequently when starting or stopping aspirin therapy.'
      }
    ],
    // Ibuprofen + Warfarin
    'Ibuprofen-Warfarin': [
      {
        type: 'Pharmacodynamic',
        severity: 'High',
        description: 'Increased risk of gastrointestinal bleeding due to combined effects on platelet function and gastric mucosa.',
        recommendation: 'Avoid concomitant use when possible. Consider alternative analgesics like acetaminophen.'
      }
    ],
    // Aspirin + Ibuprofen
    'Aspirin-Ibuprofen': [
      {
        type: 'Pharmacodynamic',
        severity: 'Moderate',
        description: 'Ibuprofen may interfere with aspirin\'s cardioprotective effects by competing for binding to COX-1.',
        recommendation: 'Take aspirin at least 30 minutes before or 8 hours after taking ibuprofen.'
      },
      {
        type: 'Adverse Effect',
        severity: 'Moderate',
        description: 'Increased risk of gastrointestinal bleeding and ulceration when NSAIDs are used together.',
        recommendation: 'Consider gastroprotective agents such as proton pump inhibitors if concurrent use is necessary.'
      }
    ],
    // Paracetamol + Warfarin
    'Paracetamol-Warfarin': [
      {
        type: 'Pharmacokinetic',
        severity: 'Low',
        description: 'Prolonged, high-dose paracetamol may enhance the effect of warfarin by affecting hepatic metabolism.',
        recommendation: 'Monitor INR more closely with prolonged or high-dose paracetamol use.'
      }
    ],
    // Caffeine + Theophylline
    'Caffeine-Theophylline': [
      {
        type: 'Pharmacokinetic',
        severity: 'Moderate',
        description: 'Competitive inhibition of metabolism leading to increased theophylline levels.',
        recommendation: 'Monitor theophylline levels and consider reducing caffeine intake.'
      },
      {
        type: 'Pharmacodynamic',
        severity: 'Moderate',
        description: 'Additive stimulant effects on central nervous system and cardiovascular system.',
        recommendation: 'Monitor for increased nervousness, insomnia, and tachycardia.'
      }
    ]
  };
  
  // Determine drug names from SMILES
  const drug1Name = knownDrugs[drug1Smiles] || 'Unknown Drug';
  const drug2Name = knownDrugs[drug2Smiles] || 'Unknown Drug';
  
  // Create a canonical key for the drug pair (alphabetical order)
  const drugPair = [drug1Name, drug2Name].sort().join('-');
  
  // Return known interactions if available
  if (knownInteractions[drugPair]) {
    return knownInteractions[drugPair];
  }
  
  // Generate random interactions for unknown pairs
  return generateRandomInteractions();
};

// Function to generate random interactions for testing
const generateRandomInteractions = (): Interaction[] => {
  const numInteractions = Math.floor(Math.random() * 3) + 1; // 1-3 random interactions
  const interactions: Interaction[] = [];
  
  const interactionTypes = ['Pharmacokinetic', 'Pharmacodynamic', 'Adverse Effect', 'Metabolic'];
  const severityLevels: ('High' | 'Moderate' | 'Low')[] = ['High', 'Moderate', 'Low'];
  
  const descriptions = [
    'May increase plasma concentrations due to CYP450 enzyme inhibition.',
    'Competitive binding to plasma proteins may increase free drug concentration.',
    'Additive effects on central nervous system depression.',
    'Increased risk of QT interval prolongation and cardiac arrhythmias.',
    'Reduced renal clearance leading to increased drug exposure.',
    'Antagonistic effects on receptor binding reducing therapeutic efficacy.',
    'Enhanced gastrointestinal adverse effects due to similar mechanisms of action.'
  ];
  
  const recommendations = [
    'Monitor drug levels and adjust dosage as needed.',
    'Consider alternative therapy with lower interaction potential.',
    'Increase monitoring frequency for adverse effects.',
    'Separate administration times to reduce interaction risk.',
    'Reduce dose of one or both medications.',
    'Consider additional prophylactic treatment to prevent adverse effects.'
  ];
  
  for (let i = 0; i < numInteractions; i++) {
    const type = interactionTypes[Math.floor(Math.random() * interactionTypes.length)];
    const severity = severityLevels[Math.floor(Math.random() * severityLevels.length)];
    const description = descriptions[Math.floor(Math.random() * descriptions.length)];
    const recommendation = Math.random() > 0.3 ? recommendations[Math.floor(Math.random() * recommendations.length)] : undefined;
    
    interactions.push({
      type,
      severity,
      description,
      recommendation
    });
  }
  
  return interactions;
};