const copyButton = document.querySelector("[data-copy]");

if (copyButton) {
  copyButton.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(copyButton.dataset.copy);
      copyButton.textContent = "Copied";
      copyButton.classList.add("copied");
      window.setTimeout(() => {
        copyButton.textContent = "Copy";
        copyButton.classList.remove("copied");
      }, 1600);
    } catch {
      copyButton.textContent = "Select command";
    }
  });
}
