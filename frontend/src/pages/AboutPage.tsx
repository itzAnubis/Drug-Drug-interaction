import React from 'react';
import { BookOpenIcon, ShieldIcon, InfoIcon, HeartIcon } from 'lucide-react';

const AboutPage = () => {
  return (
    <div className="max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">About DrugInteract</h1>
        <p className="text-gray-600">
          Learn more about our drug interaction prediction tool and how it can help healthcare professionals.
        </p>
      </div>
      
      <div className="bg-white shadow-md rounded-lg overflow-hidden mb-8">
        <div className="p-6">
          <h2 className="text-xl font-semibold text-gray-800 mb-4">Our Mission</h2>
          <p className="text-gray-600 mb-6">
            DrugInteract is designed to provide healthcare professionals and researchers with a tool to quickly identify potential interactions between medications. 
            By leveraging molecular structure data in SMILES format, we can predict interactions even for new or experimental compounds.
          </p>
          
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-8">
            <div className="bg-blue-50 p-4 rounded-lg">
              <div className="flex items-start">
                <BookOpenIcon className="h-6 w-6 text-blue-600 mt-1" />
                <div className="ml-3">
                  <h3 className="font-medium text-blue-800">Education</h3>
                  <p className="mt-1 text-sm text-gray-600">
                    Help students and professionals understand potential drug interactions based on molecular structures.
                  </p>
                </div>
              </div>
            </div>
            
            <div className="bg-teal-50 p-4 rounded-lg">
              <div className="flex items-start">
                <ShieldIcon className="h-6 w-6 text-teal-600 mt-1" />
                <div className="ml-3">
                  <h3 className="font-medium text-teal-800">Patient Safety</h3>
                  <p className="mt-1 text-sm text-gray-600">
                    Reduce adverse drug events by identifying potential interactions before they occur.
                  </p>
                </div>
              </div>
            </div>
            
            <div className="bg-amber-50 p-4 rounded-lg">
              <div className="flex items-start">
                <InfoIcon className="h-6 w-6 text-amber-600 mt-1" />
                <div className="ml-3">
                  <h3 className="font-medium text-amber-800">Research Support</h3>
                  <p className="mt-1 text-sm text-gray-600">
                    Facilitate drug discovery by predicting interactions during the development process.
                  </p>
                </div>
              </div>
            </div>
            
            <div className="bg-red-50 p-4 rounded-lg">
              <div className="flex items-start">
                <HeartIcon className="h-6 w-6 text-red-600 mt-1" />
                <div className="ml-3">
                  <h3 className="font-medium text-red-800">Clinical Decision Support</h3>
                  <p className="mt-1 text-sm text-gray-600">
                    Provide healthcare professionals with quick access to potential interaction information.
                  </p>
                </div>
              </div>
            </div>
          </div>
          
          <h2 className="text-xl font-semibold text-gray-800 mb-4">How It Works</h2>
          <p className="text-gray-600 mb-4">
            DrugInteract uses computational methods to analyze the molecular structures of drugs provided in SMILES notation. 
            The system examines structural properties, functional groups, and known interaction patterns to predict potential interactions.
          </p>
          
          <div className="bg-gray-50 p-4 rounded-lg mb-6">
            <h3 className="font-medium text-gray-800 mb-2">What is SMILES?</h3>
            <p className="text-sm text-gray-600">
              SMILES (Simplified Molecular Input Line Entry System) is a line notation for describing chemical structures. 
              It represents molecules using ASCII strings, making it easy to share and process chemical information.
              For example, aspirin is represented as <code className="bg-gray-200 px-1 py-0.5 rounded">CC(=O)OC1=CC=CC=C1C(=O)O</code>.
            </p>
          </div>
          
          <div className="border-t pt-6">
            <h2 className="text-xl font-semibold text-gray-800 mb-4">Disclaimer</h2>
            <div className="text-gray-600 text-sm">
              <p className="mb-2">
                DrugInteract is provided for educational and research purposes only. It is not intended to replace professional medical advice, diagnosis, or treatment.
              </p>
              <p className="mb-2">
                The predictions made by the system are based on computational models and may not account for all possible variables that can affect drug interactions in a clinical setting.
              </p>
              <p>
                Always consult with a qualified healthcare provider before making any decisions regarding medication use or changes to treatment plans.
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default AboutPage;