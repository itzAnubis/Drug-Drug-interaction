import React from 'react';
import { useNavigate } from 'react-router-dom';
import { ChevronRightIcon, TrashIcon, SearchIcon } from 'lucide-react';
import { useInteractionContext } from '../context/InteractionContext';

const HistoryPage = () => {
  const { interactionHistory, loadInteraction, clearHistory } = useInteractionContext();
  const [searchTerm, setSearchTerm] = React.useState('');
  const navigate = useNavigate();
  
  const filteredHistory = interactionHistory.filter(item => 
    item.drug1.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
    item.drug2.name.toLowerCase().includes(searchTerm.toLowerCase())
  );
  
  const handleSelectInteraction = (index: number) => {
    loadInteraction(index);
    navigate('/');
  };
  
  return (
    <div className="max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-gray-900 mb-2">Interaction History</h1>
        <p className="text-gray-600">
          View your previous drug interaction queries and their results.
        </p>
      </div>
      
      <div className="bg-white shadow-md rounded-lg overflow-hidden">
        <div className="p-4 border-b flex justify-between items-center">
          <div className="relative flex-grow max-w-sm">
            <div className="absolute inset-y-0 left-0 flex items-center pl-3">
              <SearchIcon className="h-5 w-5 text-gray-400" />
            </div>
            <input
              type="text"
              placeholder="Search history..."
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              className="pl-10 pr-4 py-2 border border-gray-300 rounded-md w-full focus:outline-none focus:ring-1 focus:ring-blue-500"
            />
          </div>
          
          {interactionHistory.length > 0 && (
            <button
              onClick={clearHistory}
              className="flex items-center px-3 py-2 text-sm text-red-600 hover:bg-red-50 rounded"
            >
              <TrashIcon className="h-4 w-4 mr-1" />
              Clear History
            </button>
          )}
        </div>
        
        {filteredHistory.length === 0 ? (
          <div className="p-8 text-center">
            <p className="text-gray-500">
              {interactionHistory.length === 0 
                ? "You haven't checked any drug interactions yet."
                : "No interactions match your search."}
            </p>
            {interactionHistory.length === 0 && (
              <button
                onClick={() => navigate('/')}
                className="mt-4 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Check Interactions
              </button>
            )}
          </div>
        ) : (
          <ul className="divide-y divide-gray-200">
            {filteredHistory.map((item, index) => (
              <li 
                key={index} 
                className="hover:bg-gray-50 cursor-pointer transition-colors"
                onClick={() => handleSelectInteraction(index)}
              >
                <div className="p-4 flex justify-between items-center">
                  <div>
                    <div className="flex items-center">
                      <span className="text-gray-900 font-medium">{item.drug1.name}</span>
                      <span className="mx-2 text-gray-400">+</span>
                      <span className="text-gray-900 font-medium">{item.drug2.name}</span>
                    </div>
                    <div className="mt-1 text-sm text-gray-500">
                      {new Date(item.timestamp).toLocaleString()}
                    </div>
                    <div className="mt-1 flex items-center">
                      {item.interactions.some(i => i.severity.toLowerCase() === 'high') ? (
                        <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-red-100 text-red-800">
                          High Risk
                        </span>
                      ) : item.interactions.some(i => i.severity.toLowerCase() === 'moderate') ? (
                        <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-amber-100 text-amber-800">
                          Moderate Risk
                        </span>
                      ) : (
                        <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-green-100 text-green-800">
                          Low Risk
                        </span>
                      )}
                      <span className="ml-2 text-xs text-gray-500">
                        {item.interactions.length} interaction{item.interactions.length !== 1 ? 's' : ''}
                      </span>
                    </div>
                  </div>
                  <ChevronRightIcon className="h-5 w-5 text-gray-400" />
                </div>
              </li>
            ))}
          </ul>
        )}
      </div>
    </div>
  );
};

export default HistoryPage;