import React from 'react';
import { GithubIcon, HeartIcon } from 'lucide-react';

const Footer = () => {
  return (
    <footer className="bg-white shadow-inner py-6">
      <div className="container mx-auto px-4">
        <div className="flex flex-col md:flex-row justify-between items-center">
          <div className="text-gray-600 text-sm mb-4 md:mb-0">
            © {new Date().getFullYear()} DrugInteract. All rights reserved.
          </div>
          <div className="flex items-center space-x-4">
            <span className="text-gray-600 text-sm">Made with</span>
            <HeartIcon className="h-4 w-4 text-red-500" />
            <span className="text-gray-600 text-sm">for safer prescribing</span>
            <a 
              href="https://github.com" 
              target="_blank" 
              rel="noopener noreferrer"
              className="text-gray-600 hover:text-blue-600"
            >
              <GithubIcon className="h-5 w-5" />
            </a>
          </div>
        </div>
        <div className="mt-4 text-center text-xs text-gray-500">
          This tool is for educational purposes only. Always consult healthcare professionals for medical advice.
        </div>
      </div>
    </footer>
  );
};

export default Footer;