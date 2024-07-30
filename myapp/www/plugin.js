// myplugin.js
;(function () {
  class MyPlugin {
    install() {}
    configure(pluginManager) {
      pluginManager.jexl.addFunction('formatName', feature => {
        return `<a href="${feature.name}">${feature.name}</a>`
      })
    }
  }

  // the plugin will be included in both the main thread and web worker, so
  // install plugin to either window or self (webworker global scope)
  ;(typeof self !== 'undefined' ? self : window).JBrowsePluginMyPlugin = {
    default: MyPlugin,
  }
})()