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
      status.textContent = 'Clipboard unavailable. Select the visible instructions or open their text link and copy manually.';
    }
  });
}
