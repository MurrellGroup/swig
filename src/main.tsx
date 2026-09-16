import { createRoot } from "react-dom/client";

import "./globals.css";
import "./guide-entry.css";
import SwigApp from "./swig-app";

const root = document.getElementById("root");
if (!root) throw new Error("Swig could not find its root element.");

const guideUrl = `${import.meta.env.BASE_URL}immunology/`;

createRoot(root).render(
  <>
    <nav className="guide-entry" aria-label="Documentation">
      <a href={guideUrl} target="_blank" rel="noopener noreferrer">
        B- and T-cell immunology guide
        <span className="visually-hidden"> (opens in a new tab)</span>
      </a>
      <a href={`${guideUrl}?page=13-analysis-guide`} target="_blank" rel="noopener noreferrer">
        Find background by analysis
        <span className="visually-hidden"> (opens in a new tab)</span>
      </a>
    </nav>
    <SwigApp />
  </>,
);
