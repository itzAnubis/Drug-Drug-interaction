import React from 'react';
import { Routes, Route } from 'react-router-dom';
import Header from './components/Header';
import Footer from './components/Footer';
import HomePage from './pages/HomePage';
import HistoryPage from './pages/HistoryPage';
import AboutPage from './pages/AboutPage';
import { InteractionProvider } from './context/InteractionContext';

function App() {
  return (
    <div className="min-h-screen flex flex-col bg-gray-50">
      <InteractionProvider>
        <Header />
        <main className="flex-grow container mx-auto px-4 py-8">
          <Routes>
            <Route path="/" element={<HomePage />} />
            <Route path="/history" element={<HistoryPage />} />
            <Route path="/about" element={<AboutPage />} />
          </Routes>
        </main>
        <Footer />
      </InteractionProvider>
    </div>
  );
}

export default App;