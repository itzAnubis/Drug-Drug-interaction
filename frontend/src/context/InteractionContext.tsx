import React, { createContext, useContext, useState, useEffect } from 'react';
import { InteractionResult, Drug, Interaction } from '../types/interactions';
import { predictInteractions } from '../services/interactionService';

interface InteractionContextType {
  currentInteraction: InteractionResult | null;
  interactionHistory: InteractionResult[];
  isLoading: boolean;
  checkInteractions: (drug1Name: string, drug2Name: string) => Promise<void>;
  loadInteraction: (index: number) => void;
  clearHistory: () => void;
}

const InteractionContext = createContext<InteractionContextType | undefined>(undefined);

export const useInteractionContext = () => {
  const context = useContext(InteractionContext);
  if (context === undefined) {
    throw new Error('useInteractionContext must be used within an InteractionProvider');
  }
  return context;
};

export const InteractionProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [currentInteraction, setCurrentInteraction] = useState<InteractionResult | null>(null);
  const [interactionHistory, setInteractionHistory] = useState<InteractionResult[]>([]);
  const [isLoading, setIsLoading] = useState<boolean>(false);

  useEffect(() => {
    const savedHistory = localStorage.getItem('interactionHistory');
    if (savedHistory) {
      try {
        setInteractionHistory(JSON.parse(savedHistory));
      } catch (error) {
        console.error('Failed to parse interaction history from localStorage', error);
      }
    }
  }, []);

  useEffect(() => {
    localStorage.setItem('interactionHistory', JSON.stringify(interactionHistory));
  }, [interactionHistory]);

  const checkInteractions = async (drug1Name: string, drug2Name: string) => {
    setIsLoading(true);
    
    try {
      const result = await predictInteractions(drug1Name, drug2Name);
      
      const drug1: Drug = {
        name: drug1Name,
        smiles: result.drug1_smiles,
        img_base64: result.drug1_img_base64
      };
      
      const drug2: Drug = {
        name: drug2Name,
        smiles: result.drug2_smiles,
        img_base64: result.drug2_img_base64
      };
      
      const interactionResult: InteractionResult = {
        drug1,
        drug2,
        interactions: result.interactions,
        timestamp: Date.now()
      };
      
      setCurrentInteraction(interactionResult);
      setInteractionHistory(prev => [interactionResult, ...prev]);
      
      return interactionResult;
    } catch (error) {
      console.error('Error predicting interactions:', error);
      throw error;
    } finally {
      setIsLoading(false);
    }
  };

  const loadInteraction = (index: number) => {
    if (index >= 0 && index < interactionHistory.length) {
      setCurrentInteraction(interactionHistory[index]);
    }
  };

  const clearHistory = () => {
    setInteractionHistory([]);
    setCurrentInteraction(null);
  };

  const value = {
    currentInteraction,
    interactionHistory,
    isLoading,
    checkInteractions,
    loadInteraction,
    clearHistory
  };

  return (
    <InteractionContext.Provider value={value}>
      {children}
    </InteractionContext.Provider>
  );
};