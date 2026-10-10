const path = require('path');

module.exports = {
  entry: './src/main.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
    publicPath: '/dist/',
    clean: true,
  },
  devServer: {
    static: {
      directory: __dirname,
    },
    historyApiFallback: true,
    open: true,
  },
};
