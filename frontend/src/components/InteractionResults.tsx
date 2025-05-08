import React from 'react';
import { AlertTriangleIcon, AlertCircleIcon, InfoIcon, FileTextIcon } from 'lucide-react';
import { Interaction } from '../types/interactions';

interface InteractionResultsProps {
  interactions: Interaction[];
  isLoading: boolean;
  drug1Name: string;
  drug2Name: string;
}

const InteractionResults: React.FC<InteractionResultsProps> = ({ 
  interactions, 
  isLoading, 
  drug1Name, 
  drug2Name 
}) => {
  if (isLoading) {
    return (
      <div className="bg-white shadow-md rounded-lg p-6 animate-pulse">
        <div className="h-6 bg-gray-200 rounded w-3/4 mb-4"></div>
        <div className="space-y-3">
          {[1, 2, 3].map((i) => (
            <div key={i} className="flex space-x-4">
              <div className="h-10 w-10 bg-gray-200 rounded-full"></div>
              <div className="flex-1 space-y-2 py-1">
                <div className="h-4 bg-gray-200 rounded w-5/6"></div>
                <div className="h-4 bg-gray-200 rounded w-3/4"></div>
              </div>
            </div>
          ))}
        </div>
      </div>
    );
  }

  if (!interactions.length) {
    return null;
  }

  const getSeverityIcon = (severity: string) => {
    switch (severity.toLowerCase()) {
      case 'high':
        return <AlertCircleIcon className="h-5 w-5 text-red-500" />;
      case 'moderate':
        return <AlertTriangleIcon className="h-5 w-5 text-amber-500" />;
      default:
        return <InfoIcon className="h-5 w-5 text-blue-500" />;
    }
  };

  const getSeverityClass = (severity: string) => {
    switch (severity.toLowerCase()) {
      case 'high':
        return 'bg-red-50 border-red-200';
      case 'moderate':
        return 'bg-amber-50 border-amber-200';
      case 'low':
        return 'bg-blue-50 border-blue-200';
      default:
        return 'bg-gray-50 border-gray-200';
    }
  };

  const handleExport = () => {
    // Format the interaction data
    const exportData = `
Drug Interaction Report
========================
Date: ${new Date().toLocaleString()}

Drug 1: ${drug1Name}
Drug 2: ${drug2Name}

Interactions:
${interactions.map(int => `
Severity: ${int.severity}
Type: ${int.type}
Description: ${int.description}
-------------------------`).join('\n')}
    `.trim();

    // Create a download link
    const blob = new Blob([exportData], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `interaction-report-${new Date().toISOString().slice(0, 10)}.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <div className="bg-white shadow-md rounded-lg p-6 mt-6">
      <div className="flex justify-between items-center mb-4">
        <h2 className="text-xl font-semibold text-gray-800">Interaction Results</h2>
        <button 
          onClick={handleExport}
          className="flex items-center px-3 py-1.5 text-sm text-gray-700 bg-white border border-gray-300 rounded hover:bg-gray-50"
        >
          <FileTextIcon className="h-4 w-4 mr-1" />
          Export
        </button>
      </div>
      
      <div className="mb-4 text-sm text-gray-600">
        <p>Showing potential interactions between <span className="font-medium">{drug1Name || 'Drug 1'}</span> and <span className="font-medium">{drug2Name || 'Drug 2'}</span>.</p>
      </div>
      
      <div className="space-y-4">
        {interactions.map((interaction, index) => (
          <div 
            key={index} 
            className={`p-4 rounded-md border ${getSeverityClass(interaction.severity)}`}
          >
            <div className="flex items-start">
              <div className="flex-shrink-0 mt-1">
                {getSeverityIcon(interaction.severity)}
              </div>
              <div className="ml-3">
                <h3 className="text-sm font-medium flex items-center">
                  <span className="mr-2">{interaction.type}</span>
                  <span className={`text-xs px-2 py-0.5 rounded ${
                    interaction.severity.toLowerCase() === 'high' 
                      ? 'bg-red-100 text-red-800' 
                      : interaction.severity.toLowerCase() === 'moderate'
                        ? 'bg-amber-100 text-amber-800'
                        : 'bg-blue-100 text-blue-800'
                  }`}>
                    {interaction.severity} Severity
                  </span>
                </h3>
                <div className="mt-2 text-sm text-gray-700">
                  <p>{interaction.description}</p>
                </div>
                {interaction.recommendation && (
                  <div className="mt-2 text-sm">
                    <span className="font-medium">Recommendation: </span>
                    <span className="text-gray-700">{interaction.recommendation}</span>
                  </div>
                )}
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};

export default InteractionResults;