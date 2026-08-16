import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from 'markdown-it-footnote'
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/'
  const parts = base.split('/').filter(Boolean)
  return parts.length > 0 ? `/${parts[0]}/` : '/'
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker',
  },
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  // Documenter leaves unresolved cross-references as ./@ref. The boundary API
  // page intentionally exposes a partial reference surface, so ignore only
  // those placeholders while retaining VitePress's checks for real links.
  ignoreDeadLinks: [/\/@ref/],
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', { src: `${getBaseRepository(baseTemp.base)}versions.js` }],
    ['script', { src: `${baseTemp.base}siteinfo.js` }],
    // REPLACE_ME_DOCUMENTER_VITEPRESS_NOINDEX
  ],

  markdown: {
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin)
      md.use(footnote)
      mathjax.markdownConfig(md)
    },
    theme: {
      light: 'github-light',
      dark: 'github-dark',
    },
  },
  vite: {
    plugins: [mathjax.vitePlugin],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components'),
      },
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },
  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true,
      },
    },
    nav,
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    sidebarDrawer: 'REPLACE_ME_DOCUMENTER_VITEPRESS_SIDEBAR_DRAWER',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' },
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`,
    },
  },
})
