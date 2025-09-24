// LNB-MDT 自定义JavaScript功能

document.addEventListener('DOMContentLoaded', function() {
    
    // 添加页面加载动画
    addPageLoadAnimation();
    
    // 添加代码块复制功能
    addCopyButtons();
    
    // 添加搜索增强
    enhanceSearch();
    
    // 添加滚动到顶部按钮
    addScrollToTop();
    
    // 添加代码块语法高亮增强
    enhanceCodeBlocks();
    
    // 添加表格响应式处理
    makeTablesResponsive();
    
    // 添加链接外部图标
    addExternalLinkIcons();
    
    // 添加目录导航增强
    enhanceNavigation();
    
    // 添加主题切换功能
    addThemeToggle();
    
    // 添加打印优化
    optimizePrint();
});

// 页面加载动画
function addPageLoadAnimation() {
    const content = document.querySelector('.rst-content');
    if (content) {
        content.classList.add('fade-in');
    }
}

// 添加代码块复制按钮
function addCopyButtons() {
    const codeBlocks = document.querySelectorAll('pre');
    
    codeBlocks.forEach(function(block) {
        // 跳过已经有复制按钮的代码块
        if (block.querySelector('.copybutton')) {
            return;
        }
        
        const button = document.createElement('button');
        button.className = 'copybutton';
        button.textContent = '复制';
        button.style.position = 'absolute';
        button.style.top = '10px';
        button.style.right = '10px';
        button.style.zIndex = '10';
        
        // 设置相对定位的容器
        const container = document.createElement('div');
        container.style.position = 'relative';
        block.parentNode.insertBefore(container, block);
        container.appendChild(block);
        container.appendChild(button);
        
        button.addEventListener('click', function() {
            const text = block.textContent;
            navigator.clipboard.writeText(text).then(function() {
                button.textContent = '已复制!';
                setTimeout(function() {
                    button.textContent = '复制';
                }, 2000);
            }).catch(function(err) {
                console.error('复制失败:', err);
                button.textContent = '复制失败';
                setTimeout(function() {
                    button.textContent = '复制';
                }, 2000);
            });
        });
    });
}

// 搜索增强
function enhanceSearch() {
    const searchInput = document.querySelector('input[type="text"]');
    if (searchInput) {
        searchInput.placeholder = '搜索文档...';
        
        // 添加搜索建议
        searchInput.addEventListener('input', function() {
            const query = this.value.toLowerCase();
            if (query.length > 2) {
                highlightSearchResults(query);
            }
        });
    }
}

// 高亮搜索结果
function highlightSearchResults(query) {
    const content = document.querySelector('.rst-content');
    if (!content) return;
    
    // 移除之前的高亮
    const highlights = content.querySelectorAll('.search-highlight');
    highlights.forEach(function(highlight) {
        highlight.parentNode.replaceChild(document.createTextNode(highlight.textContent), highlight);
    });
    
    // 添加新的高亮
    const walker = document.createTreeWalker(
        content,
        NodeFilter.SHOW_TEXT,
        null,
        false
    );
    
    const textNodes = [];
    let node;
    while (node = walker.nextNode()) {
        if (node.textContent.toLowerCase().includes(query)) {
            textNodes.push(node);
        }
    }
    
    textNodes.forEach(function(textNode) {
        const text = textNode.textContent;
        const regex = new RegExp(`(${query})`, 'gi');
        const highlightedText = text.replace(regex, '<span class="search-highlight">$1</span>');
        
        if (highlightedText !== text) {
            const wrapper = document.createElement('span');
            wrapper.innerHTML = highlightedText;
            textNode.parentNode.replaceChild(wrapper, textNode);
        }
    });
}

// 添加滚动到顶部按钮
function addScrollToTop() {
    const button = document.createElement('button');
    button.innerHTML = '↑';
    button.className = 'scroll-to-top';
    button.style.cssText = `
        position: fixed;
        bottom: 20px;
        right: 20px;
        width: 50px;
        height: 50px;
        border-radius: 50%;
        background: var(--primary-color);
        color: white;
        border: none;
        cursor: pointer;
        font-size: 20px;
        z-index: 1000;
        opacity: 0;
        transition: opacity 0.3s;
        box-shadow: 0 2px 10px rgba(0,0,0,0.2);
    `;
    
    document.body.appendChild(button);
    
    // 显示/隐藏按钮
    window.addEventListener('scroll', function() {
        if (window.pageYOffset > 300) {
            button.style.opacity = '1';
        } else {
            button.style.opacity = '0';
        }
    });
    
    // 点击滚动到顶部
    button.addEventListener('click', function() {
        window.scrollTo({
            top: 0,
            behavior: 'smooth'
        });
    });
}

// 代码块增强
function enhanceCodeBlocks() {
    const codeBlocks = document.querySelectorAll('pre code');
    
    codeBlocks.forEach(function(block) {
        // 添加行号
        const lines = block.textContent.split('\n');
        if (lines.length > 10) {
            const lineNumbers = document.createElement('div');
            lineNumbers.className = 'line-numbers';
            lineNumbers.style.cssText = `
                position: absolute;
                left: 0;
                top: 0;
                padding: 20px 10px;
                background: #f8f9fa;
                border-right: 1px solid #e0e0e0;
                font-family: monospace;
                font-size: 14px;
                line-height: 1.6;
                color: #666;
                user-select: none;
            `;
            
            for (let i = 1; i <= lines.length; i++) {
                const lineNumber = document.createElement('div');
                lineNumber.textContent = i;
                lineNumbers.appendChild(lineNumber);
            }
            
            const container = block.parentNode;
            container.style.position = 'relative';
            container.style.paddingLeft = '50px';
            container.appendChild(lineNumbers);
        }
        
        // 添加语言标签
        const language = block.className.match(/language-(\w+)/);
        if (language) {
            const label = document.createElement('span');
            label.className = 'language-label';
            label.textContent = language[1].toUpperCase();
            label.style.cssText = `
                position: absolute;
                top: 5px;
                right: 5px;
                background: rgba(0,0,0,0.1);
                padding: 2px 8px;
                border-radius: 3px;
                font-size: 12px;
                color: #666;
            `;
            block.parentNode.appendChild(label);
        }
    });
}

// 表格响应式处理
function makeTablesResponsive() {
    const tables = document.querySelectorAll('table');
    
    tables.forEach(function(table) {
        const wrapper = document.createElement('div');
        wrapper.className = 'table-responsive';
        wrapper.style.cssText = `
            overflow-x: auto;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        `;
        
        table.parentNode.insertBefore(wrapper, table);
        wrapper.appendChild(table);
        
        // 添加滚动提示
        if (table.scrollWidth > table.clientWidth) {
            const hint = document.createElement('div');
            hint.textContent = '← 左右滑动查看更多 →';
            hint.style.cssText = `
                text-align: center;
                color: #666;
                font-size: 12px;
                margin-top: 5px;
            `;
            wrapper.appendChild(hint);
        }
    });
}

// 添加外部链接图标
function addExternalLinkIcons() {
    const links = document.querySelectorAll('a[href^="http"]');
    
    links.forEach(function(link) {
        if (!link.querySelector('.external-icon')) {
            const icon = document.createElement('span');
            icon.className = 'external-icon';
            icon.innerHTML = '↗';
            icon.style.cssText = `
                margin-left: 5px;
                font-size: 12px;
                color: #666;
            `;
            link.appendChild(icon);
        }
    });
}

// 导航增强
function enhanceNavigation() {
    // 添加当前页面高亮
    const currentPath = window.location.pathname;
    const navLinks = document.querySelectorAll('.wy-menu-vertical a');
    
    navLinks.forEach(function(link) {
        if (link.getAttribute('href') === currentPath) {
            link.style.background = 'rgba(52, 152, 219, 0.1)';
            link.style.borderLeft = '3px solid var(--primary-color)';
        }
    });
    
    // 添加面包屑导航
    addBreadcrumb();
}

// 添加面包屑导航
function addBreadcrumb() {
    const breadcrumb = document.createElement('nav');
    breadcrumb.className = 'breadcrumb';
    breadcrumb.style.cssText = `
        background: var(--light-color);
        padding: 10px 20px;
        border-bottom: 1px solid #e0e0e0;
        font-size: 14px;
    `;
    
    const path = window.location.pathname.split('/').filter(Boolean);
    let breadcrumbHTML = '<a href="/">首页</a>';
    
    let currentPath = '';
    path.forEach(function(segment, index) {
        currentPath += '/' + segment;
        const isLast = index === path.length - 1;
        const linkText = segment.replace('.html', '').replace('-', ' ');
        
        if (isLast) {
            breadcrumbHTML += ' > <span>' + linkText + '</span>';
        } else {
            breadcrumbHTML += ' > <a href="' + currentPath + '">' + linkText + '</a>';
        }
    });
    
    breadcrumb.innerHTML = breadcrumbHTML;
    
    const content = document.querySelector('.rst-content');
    if (content) {
        content.insertBefore(breadcrumb, content.firstChild);
    }
}

// 主题切换功能
function addThemeToggle() {
    const toggle = document.createElement('button');
    toggle.innerHTML = '🌙';
    toggle.className = 'theme-toggle';
    toggle.style.cssText = `
        position: fixed;
        top: 20px;
        right: 20px;
        width: 40px;
        height: 40px;
        border-radius: 50%;
        background: var(--primary-color);
        color: white;
        border: none;
        cursor: pointer;
        font-size: 16px;
        z-index: 1000;
        box-shadow: 0 2px 10px rgba(0,0,0,0.2);
    `;
    
    document.body.appendChild(toggle);
    
    // 检查保存的主题
    const savedTheme = localStorage.getItem('theme');
    if (savedTheme === 'dark') {
        applyDarkTheme();
        toggle.innerHTML = '☀️';
    }
    
    toggle.addEventListener('click', function() {
        if (document.body.classList.contains('dark-theme')) {
            removeDarkTheme();
            toggle.innerHTML = '🌙';
            localStorage.setItem('theme', 'light');
        } else {
            applyDarkTheme();
            toggle.innerHTML = '☀️';
            localStorage.setItem('theme', 'dark');
        }
    });
}

// 应用深色主题
function applyDarkTheme() {
    document.body.classList.add('dark-theme');
    
    const style = document.createElement('style');
    style.id = 'dark-theme';
    style.textContent = `
        .dark-theme {
            background: #1a1a1a !important;
            color: #e0e0e0 !important;
        }
        .dark-theme .rst-content {
            background: #1a1a1a !important;
            color: #e0e0e0 !important;
        }
        .dark-theme .wy-side-nav-search {
            background: #2c3e50 !important;
        }
        .dark-theme .wy-menu-vertical {
            background: #2c3e50 !important;
        }
        .dark-theme .wy-menu-vertical a {
            color: #e0e0e0 !important;
        }
        .dark-theme .highlight {
            background: #2d3748 !important;
        }
        .dark-theme table {
            background: #2d3748 !important;
        }
        .dark-theme table th {
            background: #4a5568 !important;
        }
        .dark-theme table td {
            background: #2d3748 !important;
            color: #e0e0e0 !important;
        }
    `;
    document.head.appendChild(style);
}

// 移除深色主题
function removeDarkTheme() {
    document.body.classList.remove('dark-theme');
    const style = document.getElementById('dark-theme');
    if (style) {
        style.remove();
    }
}

// 打印优化
function optimizePrint() {
    // 添加打印样式
    const printStyle = document.createElement('style');
    printStyle.textContent = `
        @media print {
            .wy-side-nav,
            .wy-nav-top,
            .scroll-to-top,
            .theme-toggle,
            .breadcrumb {
                display: none !important;
            }
            .wy-content-for-nav {
                margin: 0 !important;
                padding: 0 !important;
            }
            .rst-content {
                max-width: none !important;
                padding: 0 !important;
            }
            a {
                color: black !important;
                text-decoration: underline !important;
            }
            .highlight {
                background: #f8f9fa !important;
                border: 1px solid #e0e0e0 !important;
            }
        }
    `;
    document.head.appendChild(printStyle);
}

// 添加键盘快捷键
document.addEventListener('keydown', function(e) {
    // Ctrl/Cmd + K 快速搜索
    if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
        e.preventDefault();
        const searchInput = document.querySelector('input[type="text"]');
        if (searchInput) {
            searchInput.focus();
        }
    }
    
    // ESC 清除搜索高亮
    if (e.key === 'Escape') {
        const highlights = document.querySelectorAll('.search-highlight');
        highlights.forEach(function(highlight) {
            highlight.parentNode.replaceChild(document.createTextNode(highlight.textContent), highlight);
        });
    }
});

// 添加页面加载进度条
function addLoadingProgress() {
    const progressBar = document.createElement('div');
    progressBar.className = 'loading-progress';
    progressBar.style.cssText = `
        position: fixed;
        top: 0;
        left: 0;
        width: 0%;
        height: 3px;
        background: var(--primary-color);
        z-index: 9999;
        transition: width 0.3s;
    `;
    document.body.appendChild(progressBar);
    
    // 模拟加载进度
    let progress = 0;
    const interval = setInterval(function() {
        progress += Math.random() * 30;
        if (progress >= 100) {
            progress = 100;
            clearInterval(interval);
            setTimeout(function() {
                progressBar.style.opacity = '0';
                setTimeout(function() {
                    progressBar.remove();
                }, 300);
            }, 200);
        }
        progressBar.style.width = progress + '%';
    }, 100);
}

// 初始化加载进度条
addLoadingProgress();
