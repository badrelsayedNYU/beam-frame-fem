const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {test} = require('node:test');
const root=path.join(__dirname,'..');
const html=fs.readFileSync(path.join(root,'index.html'),'utf8');
for(const file of ['fem.js','app.js','style.css']) test(`published bundle matches src/${file}`,()=>{
  assert.ok(html.includes(fs.readFileSync(path.join(root,'src',file),'utf8')),`Rebuild index.html after changing src/${file}`);
});
