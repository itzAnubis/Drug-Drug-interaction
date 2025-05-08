import React, { useState } from 'react';
import { RefreshCwIcon } from 'lucide-react';
import { useInteractionContext } from '../context/InteractionContext';

const InputForm = () => {
  const [drug1, setDrug1] = useState('');
  const [drug2, setDrug2] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  
  const { checkInteractions } = useInteractionContext();

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    
    // Validate inputs
    if (!drug1.trim() || !drug2.trim()) {
      setError('Please enter both drug names');
      return;
    }
    
    setIsLoading(true);
    
    try {
      await checkInteractions(drug1, drug2);
    } catch (err) {
      setError('An error occurred while processing your request. Please try again.');
      console.error(err);
    } finally {
      setIsLoading(false);
    }
  };
  
  const handleClear = () => {
    setDrug1('');
    setDrug2('');
    setError('');
  };

  // Example drug names
  const examples = [
    { name: 'Aspirin' },
    { name: 'Ibuprofen' },
    { name: 'Paracetamol' },
    { name: 'Warfarin' }
  ];

  return (
    <div className="bg-white shadow-md rounded-lg p-6">
      <h2 className="text-xl font-semibold text-gray-800 mb-4">Predict Drug Interactions</h2>
      
      {error && (
        <div className="bg-red-50 border-l-4 border-red-500 text-red-700 p-4 mb-4">
          <p>{error}</p>
        </div>
      )}
      
      <form onSubmit={handleSubmit}>
        <div className="space-y-4">
          <div>
            <label htmlFor="drug1" className="block text-sm font-medium text-gray-700 mb-1">
              Drug 1
            </label>
            <input
              id="drug1"
              type="text"
              value={drug1}
              onChange={(e) => setDrug1(e.target.value)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
              placeholder="Enter drug name"
            />
          </div>
          
          <div>
            <label htmlFor="drug2" className="block text-sm font-medium text-gray-700 mb-1">
              Drug 2
            </label>
            <input
              id="drug2"
              type="text"
              value={drug2}
              onChange={(e) => setDrug2(e.target.value)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
              placeholder="Enter drug name"
            />
          </div>
          
          <div className="flex items-center justify-between pt-2">
            <button
              type="button"
              onClick={handleClear}
              className="px-4 py-2 text-sm font-medium text-gray-700 bg-white border border-gray-300 rounded-md hover:bg-gray-50"
            >
              Clear
            </button>
            
            <button
              type="submit"
              disabled={isLoading}
              className="px-4 py-2 text-sm font-medium text-white bg-blue-600 border border-transparent rounded-md shadow-sm hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed flex items-center"
            >
              {isLoading ? (
                <>
                  <RefreshCwIcon className="animate-spin -ml-1 mr-2 h-4 w-4" />
                  Processing...
                </>
              ) : (
                'Check Interactions'
              )}
            </button>
          </div>
        </div>
      </form>
      
      <div className="mt-6 border-t pt-4">
        <h3 className="text-sm font-medium text-gray-700 mb-2">Example Drug Names:</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-2">
          {examples.map((example, index) => (
            <div 
              key={index}
              className="text-xs bg-gray-50 p-2 rounded cursor-pointer hover:bg-gray-100"
              onClick={() => {
                if (!drug1) setDrug1(example.name);
                else if (!drug2) setDrug2(example.name);
              }}
            >
              <div className="font-medium">{example.name}</div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};

export default InputForm;