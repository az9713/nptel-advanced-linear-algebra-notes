window.MathJax = {
  tex: {inlineMath: [['$', '$'], ['\\(', '\\)']], displayMath: [['$$', '$$'], ['\\[', '\\]']]},
  svg: {fontCache: 'local'},
  options: {enableMenu: false, skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']}
};
// An inline box with overflow other than visible loses its text baseline, so
// formulas float above the line. Only formulas wider than their block keep it.
window.fitMath = () => document.querySelectorAll('mjx-container:not([display=true])').forEach(box => {
  box.classList.remove('fits');
  const svg = box.querySelector('svg');
  const block = box.closest('p,li,td,th,h1,h2,h3,figcaption,summary,div');
  const width = svg ? svg.getBoundingClientRect().width : 0;
  if (width && block && width <= block.clientWidth) box.classList.add('fits');
});
window.MathJax.startup = {pageReady: () => MathJax.startup.defaultPageReady().then(window.fitMath)};
addEventListener('resize', window.fitMath);
addEventListener('toggle', window.fitMath, true);
