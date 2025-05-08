import React, { useState } from 'react';
import InputForm from '../components/InputForm';
import MoleculeViewer from '../components/MoleculeViewer';
import InteractionResults from '../components/InteractionResults';
import { useInteractionContext } from '../context/InteractionContext';

const HomePage = () => {
  const { currentInteraction, isLoading } = useInteractionContext();
  
  return (
    <div className="max-w-6xl mx-auto">
      <div className="text-center mb-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">Drug Interaction Predictor</h1>
        <p className="text-gray-600 max-w-3xl mx-auto">
          Enter the SMILES notation for two drugs to predict potential interactions between them.
          Our model analyzes molecular structures to identify possible adverse effects.
        </p>
      </div>
      
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        <div className="lg:col-span-1">
          <InputForm />
        </div>
        
        <div className="lg:col-span-2">
          {currentInteraction && (
            <>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
                <MoleculeViewer 
                  smiles={currentInteraction.drug1.smiles} 
                  name={currentInteraction.drug1.name || 'Drug 1'} 
                />
                <MoleculeViewer 
                  smiles={currentInteraction.drug2.smiles} 
                  name={currentInteraction.drug2.name || 'Drug 2'} 
                />
              </div>
              
              <InteractionResults 
                interactions={currentInteraction.interactions} 
                isLoading={isLoading}
                drug1Name={currentInteraction.drug1.name || 'Drug 1'}
                drug2Name={currentInteraction.drug2.name || 'Drug 2'}
              />
            </>
          )}
          
          {!currentInteraction && !isLoading && (
            <div className="bg-white shadow-md rounded-lg p-6 h-full flex flex-col items-center justify-center text-center">
              <div className="p-6 max-w-lg">
                <h3 className="text-xl font-medium text-gray-900 mb-2">How It Works</h3>
                <p className="text-gray-600 mb-4">
                  Our tool uses advanced algorithms to predict potential interactions between drugs based on their molecular structures.
                </p>
                <ol className="text-left text-gray-600 space-y-3 mb-6">
                  <li className="flex">
                    <span className="flex-shrink-0 w-6 h-6 bg-blue-100 text-blue-800 rounded-full flex items-center justify-center mr-2">1</span>
                    <span>Enter the SMILES notation for two drugs</span>
                  </li>
                  <li className="flex">
                    <span className="flex-shrink-0 w-6 h-6 bg-blue-100 text-blue-800 rounded-full flex items-center justify-center mr-2">2</span>
                    <span>Our model analyzes molecular structures</span>
                  </li>
                  <li className="flex">
                    <span className="flex-shrink-0 w-6 h-6 bg-blue-100 text-blue-800 rounded-full flex items-center justify-center mr-2">3</span>
                    <span>Review predicted interactions and severity levels</span>
                  </li>
                </ol>
                <p className="text-sm text-gray-500">
                  Not sure where to find SMILES notations? Use the example drugs provided in the input form.
                </p>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default HomePage;