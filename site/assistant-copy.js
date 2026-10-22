'use strict';
for (const button of document.querySelectorAll('[data-copy-instructions]')) {
  button.addEventListener('click', async () => {
    const region = document.getElementById(button.dataset.copyInstructions);
    const code = region.querySelector('code') || region;
    const status = document.querySelector('[data-copy-status]');
    try {
      await navigator.clipboard.writeText(code.textContent.trim() + '\n');
      status.textContent = 'Complete assistant instructions copied. Paste them into your own conversation.';
    } catch (error) {
      const details = region.closest?.('details');
      if (details) details.open = true;
      status.textContent = 'Clipboard unavailable. Use the expanded full prompt or download the text and copy manually.';
    }
  });
}
