import { createRoot } from "react-dom/client";

import "./globals.css";
import SwigApp from "./swig-app";

const root = document.getElementById("root");
if (!root) throw new Error("Swig could not find its root element.");

createRoot(root).render(<SwigApp />);

