import React from 'react';

interface MoleculeViewerProps {
  smiles: string;
  name: string;
  img_base64?: string;
}

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ smiles, name, img_base64 }) => {
  return (
    <div className="flex flex-col items-center p-4 bg-white rounded-lg shadow">
      <h3 className="text-md font-medium text-gray-800 mb-2">{name}</h3>
      <div className="relative w-full aspect-square bg-gray-100 rounded flex items-center justify-center overflow-hidden">
        {img_base64 ? (
          <img 
            src={`data:image/png;base64,${img_base64}`}
            alt={`Molecular structure of ${name}`}
            className="max-w-full max-h-full object-contain"
          />
        ) : (
          <img 
            src={`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(smiles)}/image`}
            alt={`Molecular structure of ${name}`}
            className="max-w-full max-h-full object-contain"
            onError={(e) => {
              const target = e.target as HTMLImageElement;
              target.onerror = null;
              target.src = 'https://placehold.co/400x400?text=Molecule+Visualization+Unavailable';
            }}
          />
        )}
      </div>
      <div className="mt-2 w-full">
        <div className="text-xs text-gray-500 break-all">
          <span className="font-medium">SMILES:</span> {smiles}
        </div>
      </div>
    </div>
  );
};

export default MoleculeViewer;