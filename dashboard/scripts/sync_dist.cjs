const fs = require("fs");
const path = require("path");
const ts = require("typescript");

const rootDir = path.resolve(__dirname, "..");
const srcDir = path.join(rootDir, "src");
const distDir = path.join(rootDir, "dist");
const assetsDir = path.join(distDir, "assets");

function readSourceMain() {
  let source = fs.readFileSync(path.join(srcDir, "main.ts"), "utf8");
  source = source.replace(/^import\s+"\.\/styles\.css";\s*/m, "");
  source = source.replace(/declare global \{[\s\S]*?\n\}\n\n/m, "");
  source = source.replace(/\nexport \{\};?\s*$/m, "\n");
  return source;
}

function buildAppJs() {
  const transpiled = ts.transpileModule(readSourceMain(), {
    compilerOptions: {
      target: ts.ScriptTarget.ES2020,
      module: ts.ModuleKind.None,
      removeComments: false,
    },
  }).outputText.trim();

  return `(function () {\n${transpiled}\n})();\n`;
}

function buildIndexHtml() {
  return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>TASOPT Diagnostic Dashboard</title>
  <link rel="stylesheet" href="./assets/app.css" />
  <script>
    window.MathJax = {
      tex: {
        inlineMath: [["\\\\(", "\\\\)"], ["$", "$"]],
      },
      chtml: {
        scale: 0.96,
      },
      startup: {
        ready() {
          MathJax.startup.defaultReady();
          document.dispatchEvent(new Event("mathjax-ready"));
        },
      },
    };
    window.PlotlyConfig = {
      MathJaxConfig: "local",
    };
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
  <script defer src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
</head>
<body>
  <div id="app"></div>
  <script src="./assets/app.js"></script>
</body>
</html>
`;
}

function main() {
  fs.rmSync(distDir, { recursive: true, force: true });
  fs.mkdirSync(assetsDir, { recursive: true });
  fs.writeFileSync(path.join(assetsDir, "app.js"), buildAppJs());
  fs.copyFileSync(path.join(srcDir, "styles.css"), path.join(assetsDir, "app.css"));
  fs.writeFileSync(path.join(distDir, "index.html"), buildIndexHtml());
}

main();
