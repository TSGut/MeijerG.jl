document.addEventListener("DOMContentLoaded", () => {
    const html = document.documentElement;

    // 1. Remove all theme classes
    html.classList.remove(
        "theme--documenter-light",
        "theme--documenter-dark",
        "theme--catppuccin-latte",
        "theme--catppuccin-frappe",
        "theme--catppuccin-macchiato",
        "theme--catppuccin-mocha"
    );

    // 2. Set the default theme = Catppuccin Macchiato
    html.classList.add("theme--catppuccin-macchiato");

    // 3. Modify the theme switcher to include only:
    //    - Macchiato (default)
    //    - Latte (optional)
    const menu = document.querySelector(".theme-selector");
    if (menu) {
        // Clear all existing entries
        while (menu.firstChild) menu.removeChild(menu.firstChild);

        // Create Latte option
        const latte = document.createElement("button");
        latte.textContent = "Catppuccin Latte";
        latte.onclick = () => {
            html.classList.remove("theme--catppuccin-macchiato");
            html.classList.add("theme--catppuccin-latte");
            localStorage.setItem("documenter-theme", "catppuccin-latte");
        };

        // Create Macchiato option
        const macchiato = document.createElement("button");
        macchiato.textContent = "Catppuccin Macchiato";
        macchiato.onclick = () => {
            html.classList.remove("theme--catppuccin-latte");
            html.classList.add("theme--catppuccin-macchiato");
            localStorage.setItem("documenter-theme", "catppuccin-macchiato");
        };

        menu.appendChild(macchiato);
        menu.appendChild(latte);
    }

    // 4. Enforce default on first load unless user explicitly switched
    const stored = localStorage.getItem("documenter-theme");
    if (!stored) {
        // no manual user choice → enforce default
        localStorage.setItem("documenter-theme", "catppuccin-macchiato");
    } else {
        // user previously chose Latte → respect that
        html.classList.remove("theme--catppuccin-macchiato", "theme--catppuccin-latte");
        html.classList.add(`theme--${stored}`);
    }
});
